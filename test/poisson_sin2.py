import test_tools
import fempy
import fempy.mesh
import fempy.geometry
import fempy.expression_operators as eo
import pylinear.array as num



sol = (eo.SIN, (eo.TIMES, 5, (eo.TIMES, (eo.POWER,(eo.VARIABLE,"0"),2), 
                              (eo.POWER,(eo.VARIABLE,"1"),2))))
grid_vectors = [num.array([2,0], num.Float),
                num.array([0,2], num.Float)]
mesh = fempy.mesh.tTwoDimensionalMesh(
  [fempy.mesh.tShapeSection(fempy.geometry.getParallelogram(grid_vectors), "dirichlet")])
test_tools.adaptiveDemo(sol, mesh, max_iterations = 20)
