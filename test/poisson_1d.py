import math

import pylinear.array as num

import fempy
import fempy.mesh
import fempy.geometry
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

mesh = fempy.mesh.OneDimensionalMesh(-1, 1, 200, "dirichlet", "dirichlet")
constraints = fempy.solver.get_dirichlet_constraints(mesh, u_d)
solution = fempy.solver.solve_poisson(mesh, f, constraints)
fempy.visualization.visualize_1d_mesh_function(solution, ",,1d_solution.data")
