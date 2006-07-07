import pymbolic
import test_tools
import fempy
import fempy.mesh
import fempy.geometry
import pylinear.array as num

sol = pymbolic.parse("math.sin(5*x[0]**2*x[1]**2)")
grid_vectors = [num.array([2,0], num.Float),
                num.array([0,2], num.Float)]
mesh = fempy.mesh.TwoDimensionalMesh(
  [fempy.mesh.ShapeSection(fempy.geometry.get_parallelogram(grid_vectors), "dirichlet")])
test_tools.adaptive_demo(sol, mesh, max_iterations=20)
