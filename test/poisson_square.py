import math

import pylinear.matrices as num

import fempy
import fempy.mesh
import fempy.geometry
import fempy.expression as expression

import test_tools

case = "zero_boundary"

if case == "zero_boundary":
  expr_x = ("-", 1, ("**",("variable","0"), 2))
  expr_y = ("-", 1, ("**",("variable","1"), 2))
  sol = ("*",expr_x, expr_y)
elif case == "constant_rhs":
  sol = ("+", ("**", ("variable","0"), 2), ("**", ("variable", "1"), 2)) 

grid_vectors = [num.array([2,0], num.Float),
                num.array([0,2], num.Float)]
mesh = fempy.mesh.tTwoDimensionalMesh(
  [fempy.mesh.tShapeSection(fempy.geometry.getParallelogram(grid_vectors), "dirichlet")])
test_tools.adaptiveDemo(sol, mesh, max_iterations = 10)
