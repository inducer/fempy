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

def f(x):
  if x[0] > 0:
    return -500
  else:
    return 500

def u_d(x):
  if x[0] > 0:
    return 0
  else:
    return 0

mesh = fempy.mesh.tOneDimensionalMesh(-1, 1, 200, "dirichlet", "dirichlet")
constraints = fempy.solver.getDirichletConstraints(mesh, u_d)
solution = fempy.solver.solvePoisson(mesh, f, constraints)
fempy.visualization.visualize1DMeshFunction(solution, ",,1d_solution.data")
