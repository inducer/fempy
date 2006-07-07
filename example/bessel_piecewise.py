import pylinear.operation as op
import fempy.mesh
import fempy.geometry as geometry
import fempy.solver as solver
import fempy.visualization as visualization

def needs_refinement(vert_origin, vert_destination, vert_apex, area):
    return area >= 3e-2

alpha_1 = 1; alpha_2 = 40

def alpha(x):
    if op.norm_2(x) < 0.5:
        return alpha_1
    else:
        return alpha_2

mesh = fempy.mesh.TwoDimensionalMesh(
    [fempy.mesh.ShapeSection(fempy.geometry.get_circle(1), "dirichlet"),
     fempy.mesh.ShapeSection(fempy.geometry.get_circle(0.5), "unconstrained")], 
    refinement_func = needs_refinement)
constraints = solver.get_dirichlet_constraints(mesh, u_d = lambda x: 0)
eigensolver = solver.LaplacianEigenproblemSolver(
  mesh, constrained_nodes = constraints, g = alpha)
eigensolver.setup_constraints(constraints)
solutions = eigensolver.solve(0, tolerance = 1e-10, number_of_eigenvalues = 20)

for evalue, emode in solutions:
    visualization.visualize("vtk", (",,result.vtk", ",,result_grid.vtk"), emode.real)
    if raw_input("[enter for next, q for quit]") == "q":
        break
