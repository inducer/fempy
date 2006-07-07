import sys, math
import pylinear.array as num
import pylinear.operation as op
import fempy.mesh
import fempy.geometry as geometry
import fempy.solver as solver
import fempy.visualization as visualization
import fempy.element_norm as element_norm
import fempy.mesh_function as mesh_function
import fempy.tools as tools
import scipy.special as special
import fempy.eoc as eoc
import scipy.optimize

def needsRefinement( vert_origin, vert_destination, vert_apex, area ):
    return area >= max_tri_area

alpha_1 = 40; alpha_2 = 1
origin = num.zeros((2,), num.Float)

def alpha(r):
    if op.norm_2(r) < 0.5:
        return alpha_1
    else:
        return alpha_2

def R(lambd, n, c11, c21, c22, r):
    if r < 0.5:
        return R_1(lambd, n, c11, r)
    else:
        return R_2(lambd, n, c21, c22, r)

def R_1(lambd, n, c11, r):
    rho = math.sqrt(lambd*alpha_1) * r
    return c11 * special.jn(n, rho)

def R_1diff(lambd, n, c11, r):
    rho = math.sqrt(lambd*alpha_1) * r
    return c11 * special.jvp(n, rho)\
           * math.sqrt(lambd*alpha_1)

def R_2(lambd, n, c21, c22, r):
    rho = math.sqrt(lambd*alpha_2) * r
    return c21 * special.jn(n, rho) + c22 * special.yn(n, rho)

def R_2diff(lambd, n, c21, c22, r):
    rho = math.sqrt(lambd*alpha_2) * r
    return (c21 * special.jvp(n, rho) + c22 * special.yvp(n, rho)) \
           * math.sqrt(lambd*alpha_2)

# Are we using the right functions?
#tools.write1DGnuplotGraphs(
    #lambda x: (special.yn(1, x), special.yvp(1, x)),
    #a = 0.5, b = 30)
#sys.exit()

do_visualization = raw_input("visualize? [n]") == "y"

eigenvalue_eoc = eoc.EOCRecorder()
eigenfunc_eoc = eoc.EOCRecorder()
max_tri_area = 8e-2
for step in range(4):
    max_tri_area *= 0.8

    mesh = fempy.mesh.tTwoDimensionalMesh(
        [fempy.mesh.tShapeSection(fempy.geometry.getCircle(1), "dirichlet"),
         fempy.mesh.tShapeSection(fempy.geometry.getCircle(0.5), "unconstrained")], 
        refinement_func = needsRefinement)
    print "ITERATION %d: %d elements, %d nodes" % (step, len(mesh.elements()),
                                                   len(mesh.dofManager()))
    constraints = solver.getDirichletConstraints(mesh, u_d = lambda x: 0)
    eigensolver = solver.tLaplacianEigenproblemSolver(
      mesh, constrained_nodes = constraints, g = alpha)
    eigensolver.setupConstraints(constraints)
    solutions = eigensolver.solve(0, tolerance = 1e-10, number_of_eigenvalues = 20)

    evalue_error = 0
    for n, (fempy_evalue, fempy_emode) in enumerate(solutions[0:1]):
        if fempy_evalue.imag != 0:
            continue

        def findMyZero(arry):
            mu, nu, lambd = arry
            lambd /= 100
            c11 = 1
            c21 = mu
            c22 = nu
            #print "trying ", mu, nu, lambd
            return [R_1(lambd, n, c11, 0.5)-R_2(lambd, n, c21, c22, 0.5),
                   R_1diff(lambd, n, c11, 0.5) - R_2diff(lambd, n, c21, c22, 0.5),
                   R_2(lambd, n, c21, c22, 1)]
        mu_0, nu_0, lambda_0 = scipy.optimize.fsolve(
            findMyZero, [0, 0, fempy_evalue.real * 100])
        c11_0 = 1
        c21_0 = mu_0
        c22_0 = nu_0
        lambda_0 /= 100
        evalue_error += abs(lambda_0-fempy_evalue)**2
        print "found", mu_0, nu_0, lambda_0, fempy_evalue
        print "f residual:", op.norm_2(num.array(findMyZero([mu_0, nu_0, lambda_0*100])))

        def Rvector(r_vec):
            return R(lambda_0, n, c11_0, c21_0, c22_0, op.norm_2(r_vec))
 
        if do_visualization:
            tools.write_1d_gnuplot_graph(
                lambda r: R(lambda_0, n, c11_0, c21_0, c22_0, r), 
                a = 0.01, b = 1)
            visualization.visualize("vtk", (",,result.vtk", ",,result_grid.vtk"), fempy_emode.real)
            raw_input("[showing computed, enter]")

            ana = mesh_function.discretizeFunction(mesh, Rvector)
            visualization.visualize("vtk", (",,result.vtk", ",,result_grid.vtk"), ana)
            raw_input("[showing analytic, enter]")

        fempy_emode *=  Rvector(origin) / fempy_emode(origin)

        my_estimator = element_norm.makeL2ErrorNormSquared(Rvector, fempy_emode)
        total_error = tools.sum_over(my_estimator, mesh.elements())

        eigenfunc_eoc.add_data_point(math.sqrt(len(mesh.elements())),
                                     math.sqrt(abs(total_error)))

    eigenvalue_eoc.add_data_point(math.sqrt(len(mesh.elements())),
                                  math.sqrt(evalue_error))

print "-------------------------------------------------------"
print "Eigenvalue EOC overall:", eigenvalue_eoc.estimate_order_of_convergence()[0,1]
print "EOC Gliding means:"
gliding_means = eigenvalue_eoc.estimate_order_of_convergence(3)
gliding_means_iterations,dummy = gliding_means.shape
for i in range(gliding_means_iterations):
    print "Iteration %d: %f" % (i, gliding_means[i,1])
print "-------------------------------------------------------"
eigenvalue_eoc.write_gnuplot_file(",,eigenvalue_conv.data")
print "-------------------------------------------------------"
print "Eigenfunction EOC overall:", eigenfunc_eoc.estimate_order_of_convergence()[0,1]
print "EOC Gliding means:"
gliding_means = eigenfunc_eoc.estimate_order_of_convergence(3)
gliding_means_iterations,dummy = gliding_means.shape
for i in range(gliding_means_iterations):
    print "Iteration %d: %f" % (i, gliding_means[i,1])
print "-------------------------------------------------------"
eigenfunc_eoc.write_gnuplot_file(",,eigenfunc_conv.data")
