import math

import pylinear.matrices as num

import fempy
import fempy.mesh
import fempy.geometry
import fempy.expression as expression
import fempy.expression_operators as eo
import fempy.visualization
import fempy.solver

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

mesh = fempy.mesh.tOneDimensionalMesh(-1, 1, 1000, "dirichlet", "dirichlet")

def f(x):
  if abs(x[0]) < 0.005:
    return 0
  else:
    return 500 * math.sin(6/x[0])

def u_d(x):
  if x[0] > 0:
    return 3
  else:
    return -3

constraints = fempy.solver.getDirichletConstraints(mesh, u_d)
solution = fempy.solver.solvePoisson(mesh, f,
                                     constraints)

fempy.visualization.visualize1DMeshFunction(solution, ",,1d_solution.data")
