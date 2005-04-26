import math

import pylinear.array as num

import fempy
import fempy.mesh
import fempy.geometry
import fempy.expression as expression
import fempy.expression_operators as eo

import test_tools

#import psyco

#psyco.log()
#psyco.profile()

case = "zero_boundary"

if case == "zero_boundary":
  expr_x = (eo.MINUS, 1, (eo.POWER,(eo.VARIABLE,"0"), 2))
  expr_y = (eo.MINUS, 1, (eo.POWER,(eo.VARIABLE,"1"), 2))
  sol = (eo.TIMES,expr_x, expr_y)
elif case == "constant_rhs":
  sol = (eo.PLUS, (eo.POWER, (eo.VARIABLE,"0"), 2), (eo.POWER, (eo.VARIABLE, "1"), 2)) 

grid_vectors = [num.array([2,0], num.Float),
                num.array([0,2], num.Float)]
mesh = fempy.mesh.tTwoDimensionalMesh(
  [fempy.mesh.tShapeSection(fempy.geometry.getParallelogram(grid_vectors), "dirichlet")])
test_tools.adaptiveDemo(sol, mesh, max_iterations = 13)
