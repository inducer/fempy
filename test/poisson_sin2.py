import test_tools
import fempy
import fempy.mesh
import fempy.geometry



sol = ("sin", ("*", 5, ("*", ("**",("variable","0"),2), ("**",("variable","1"),2))))
mesh = fempy.mesh.tTwoDimensionalMesh(
  fempy.geometry.getUnitCellGeometry(edge_length = 2))
test_tools.adaptiveDemo(sol, mesh, max_iterations = 3)
