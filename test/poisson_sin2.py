import test_tools
import fempy
import fempy.mesh
import fempy.geometry
import pylinear.matrices as num



sol = ("sin", ("*", 5, ("*", ("**",("variable","0"),2), ("**",("variable","1"),2))))
grid_vectors = [num.array([2,0], num.Float),
                num.array([0,2], num.Float)]
mesh = fempy.mesh.tTwoDimensionalMesh(
  fempy.geometry.getUnitCellGeometry(grid_vectors, 0.6))
test_tools.adaptiveDemo(sol, mesh, max_iterations = 3)
