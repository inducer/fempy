import math

import test_tools
import fempy
import fempy.mesh
import fempy.geometry
import fempy.expression as expression

case = "zero_boundary"

if case == "zero_boundary":
  expr_x = ("-", 1, ("**",("variable","0"), 2))
  expr_y = ("-", 1, ("**",("variable","1"), 2))
  sol = ("*",expr_x, expr_y)
elif case == "constant_rhs":
  sol = ("+", ("**", ("variable","0"), 2), ("**", ("variable", "1"), 2)) 

mesh = fempy.mesh.tTwoDimensionalMesh(
  [fempy.mesh.tShapeSection(fempy.geometry.getParallelogram(edge_length = 2), True)])
test_tools.adaptiveDemo(sol, mesh, max_iterations = 6)
