#! /usr/bin/python

import profile
import pstats

# System imports --------------------------------------------------------------
import math
import sys

# Numerics imports ------------------------------------------------------------
import pylinear.matrices as num
import pylinear.linear_algebra as la

# FEM imports -----------------------------------------------------------------
import pyangle
import visualization
import spatial_btree
from matrix_builder import *
from element import *
from stopwatch import *
from solver import *
from tools import *
from geometry import *
import eoc




# driver ----------------------------------------------------------------------
def makeSolutionFunction(elements, solution, finder = None):
  if finder is None:
    finder = spatial_btree.buildElementFinder(elements)
  def f(point):
    el = finder(point)
    if el is None:
      raise RuntimeError, "Solution not defined at %s" % str(point)
    return el.getSolutionFunction(solution)(point)
  return f




def solveAdaptively(mesh, solve, mark_for_refinement, iteration_limit = None):
  mesh.generate()

  iteration_count = 0
  solution = None
  while iteration_count < iteration_limit:
    refinement_info, solution = solve(mesh, solution)
    mesh_change = mesh.getRefinement(lambda el: mark_for_refinement(refinement_info, el))
    if mesh_change is None:
      return solution

    mesh = mesh_change.meshAfter()
    mesh.generate()
    solution = mesh_change.changeSolution(solution)

    iteration_count += 1



  
def junk():
  def needsRefinement(vert_origin, vert_destination, vert_apex, area):
    return area > 0.001
  shape = [ (0.01,0), (1,0), (1,1), (0,1) ]
  nodes, elements = buildShapeGeometry(dof_manager, shape, needsRefinement, False)




def poissonDemo(n):
  center = num.array([0.5, 0.5])

  def f(x):
    if norm2(x - center) < 0.3:
      return -20
    else:
      return 20

  def u_d(x):
    if 0.1 < x[0] < 0.9 and x[1] < 0.5:
      return 1
    else:
     return 0

  job = tJob("geometry")
  mesh = tRectangularMesh(n, n)
  nodes, elements, hanging_nodes = mesh.getData()
  job.done()

  # make the edge nodes dirichlet nodes
  def isEdgeNode(node):
    x = node.coordinates()
    return x[0] in [0,1] or x[1] in [0,1]

  dirichlet_nodes = filter(isEdgeNode, nodes)
  
  solution = solvePoisson(mesh.dofManager(), elements, dirichlet_nodes, f, u_d)

  #s_f1 = makeSolutionFunction(elements, solution, finder)

  visualization.writeMatlabFile("/tmp/visualize.m", mesh.dofManager(), elements, solution)
  #visualization.writeGnuplotFile("+result.dat", dof_manager, elements, solution)
  #visualization.writeVtkFile("+result.vtk", dof_manager, elements, solution)
  



def poissonTest():
  width = 1.
  height = 1.
  
  nx = 40
  ny = 40

  a = 5

  def f(point):
    x = point[0]
    y = point[1]
    return \
      -4 * a **2 * (x**2 * y**4 + x**4 * y**2) * math.sin(a * x**2 * y**2) + \
      2 * a * (y**2 + x**2) * math.cos(a * x**2 * y**2)

  def solution(point):
    x = point[0]
    y = point[1]
    return math.sin(a * x**2 * y** 2)
    
  dof_manager = tDOFManager()
  nodes, elements = buildRectangularGeometry(dof_manager, width / nx, height / ny, nx, ny, True)

  # make the edge nodes dirichlet nodes
  def isEdgeNode(node):
    x = node.coordinates()
    return x[0] in [0,1] or x[1] in [0,1]

  dirichlet_nodes = filter(isEdgeNode, nodes)
  solution_vector = solvePoisson(dof_manager, elements, dirichlet_nodes, f, solution)

  visualization.writeVtkFile("+result.vtk", dof_manager, elements, solution_vector)

  def errorFunctionL2(point, solution_func_value):
    result = (solution(point) - solution_func_value) ** 2
    assert result > 0 
    return result

  error_integral = 0
  for el in elements:
    node_values = num.take(solution_vector, el.nodeNumbers())
    error_integral += el.getVolumeIntegralOver(errorFunctionL2, node_values)

  print "L_2 error estimate:", error_integral
    


#poissonTest()
poissonDemo(10)

def doProfile():
  profile.run("poissonDemo()", ",,profile")
  p = pstats.Stats(",,profile")
  p.sort_stats("time").print_stats()


