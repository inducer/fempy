import pylinear.matrix_tools as mtools
import fempy.mesh as mesh
import fempy.geometry as geometry
import fempy.solver as solver
import fempy.visualization as visualization

def needsRefinement( vert_origin, vert_destination, vert_apex, area ):
    return area >= 3e-2

def rhscoeff(x):
    if mtools.norm2(x) < 0.5:
        return 40
    else:
        return 1

mesh = fempy.mesh.tTwoDimensionalMesh(
    [mesh.tShapeSection(fempy.geometry.getCircle(1), "dirichlet"),
     mesh.tShapeSection(fempy.geometry.getCircle(0.5), "unconstrained")], 
    refinement_func = needsRefinement)
constraints = solver.getDirichletConstraints(mesh, u_d = lambda x: 0)
eigensolver = solver.tLaplacianEigenproblemSolver(
  mesh, constrained_nodes = constraints, g = rhscoeff)
eigensolver.setupConstraints(constraints)
solutions = eigensolver.solve(0, tolerance = 1e-10, number_of_eigenvalues = 20)

for evalue, emode in solutions:
    visualization.visualize("vtk", (",,result.vtk", ",,result_grid.vtk"), emode.real)
    if raw_input("[enter for next, q for quit]") == "q":
        break
