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
import fempy.visualization as vis
import fempy.eoc as eoc
import fempy.error_estimator as error_estimator
from fempy.matrix_builder import *
from fempy.element import *
from fempy.stopwatch import *
from fempy.solver import *
from fempy.tools import *
from fempy.mesh import *




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




def solveAdaptively(mesh, solve, element_needs_refining, iteration_limit = None):
  mesh.generate()

  iteration_count = 0
  solution = None
  while True:
    refinement_info, solution = solve(mesh, solution)

    iteration_count += 1
    if iteration_limit is not None and iteration_count >= iteration_limit:
      return mesh, solution

    mesh_change = mesh.getRefinement(lambda el: element_needs_refining(refinement_info, el))

    if mesh_change is None:
      return mesh,solution

    mesh = mesh_change.meshAfter()
    mesh.generate()
    solution = mesh_change.changeSolution(solution)



  
def adaptiveDemo():
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

  job = tJob("geometry")
  mesh = tShapedMesh([(0,0), (0,1), (1,1), (1,0)])
  mesh.generate()
  job.done()

  def isEdgeNode(node):
    x = node.coordinates()
    return x[0] in [0,1] or x[1] in [0,1]

  def solve(new_mesh, start_solution_vector):
    dirichlet_nodes = filter(isEdgeNode, new_mesh.nodes())
    solution_vector = solvePoisson(new_mesh, dirichlet_nodes, f, solution, start_solution_vector)

    my_estimator = error_estimator.tAnalyticSolutionErrorEstimator(new_mesh, solution_vector,
      solution)

    max_error = max(map(my_estimator, new_mesh.elements()))
    return (my_estimator, max_error), solution_vector

  def decideOnRefinement((error_estimator, max_error), element):
    return error_estimator(element) > 0.5 * max_error

  new_mesh, solution_vector = solveAdaptively(mesh, solve, decideOnRefinement, 15)

  vis.writeMatlabFile("/tmp/visualize.m", new_mesh.dofManager(), new_mesh.elements(), solution_vector)
  #vis.writeGnuplotFile("+result.dat", dof_manager, elements, solution)
  #vis.writeVtkFile("+result.vtk", dof_manager, elements, solution)




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
  mesh.generate()
  nodes = mesh.nodes()
  job.done()

  # make the edge nodes dirichlet nodes
  
  solution = solvePoisson(mesh.dofManager(), mesh.elements(), dirichlet_nodes, f, u_d)

  #s_f1 = makeSolutionFunction(elements, solution, finder)

  vis.writeMatlabFile("/tmp/visualize.m", mesh.dofManager(), mesh.elements(), solution)
  



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

  vis.writeVtkFile("+result.vtk", dof_manager, elements, solution_vector)

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
#poissonDemo(20)
adaptiveDemo()

def doProfile():
  profile.run("poissonDemo()", ",,profile")
  p = pstats.Stats(",,profile")
  p.sort_stats("time").print_stats()


