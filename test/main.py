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
  iteration_count = 0
  solution = None
  while True:
    print "--------------------------------------------------------------------"
    print "ITERATION %d: %d elements, %d nodes" % (iteration_count, len(mesh.elements()), len(mesh.nodes()))
    print "--------------------------------------------------------------------"
    refinement_info, solution = solve(mesh, solution)

    iteration_count += 1
    if iteration_limit is not None and iteration_count >= iteration_limit:
      return mesh, solution

    job = tJob("refining...")
    mesh_change = mesh.getRefinement(lambda el: element_needs_refining(refinement_info, el))

    if mesh_change is None:
      job.done()
      return mesh,solution

    mesh = mesh_change.meshAfter()
    mesh.generate()
    job.done()

    #job = tJob("adapting...")
    #solution = mesh_change.changeSolution(solution)
    #job.done()

    solution = None



  
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
  #boundary = [(0,0), (0,1), (1,1), (1,0)] # square
  boundary = [(0,0), (1,0), (2,1), (1,1)] # parallelogram
  
  mesh = tTwoDimensionalShapedMesh(boundary)
  mesh.generate()
  job.done()

  eoc_rec = eoc.tEOCRecorder()

  def solve(new_mesh, start_solution_vector = None):
    if start_solution_vector is None:
      start_solution_vector = num.zeros((new_mesh.dofManager().countDegreesOfFreedom(),), num.Float)

    solution_vector = solvePoisson(new_mesh, new_mesh.boundaryNodes(),
      f, solution, start_solution_vector)

    my_estimator = error_estimator.tAnalyticSolutionL2ErrorEstimator(
      new_mesh, solution_vector, solution)

    max_error = max(map(my_estimator, new_mesh.elements()))
    eoc_rec.addDataPoint(len(new_mesh.elements())**0.5, my_estimator.estimateTotalError())
    return (my_estimator, max_error), solution_vector

  def decideOnRefinement((error_estimator, max_error), element):
    return error_estimator(element) > 0.5 * max_error

  new_mesh, solution_vector = solveAdaptively(mesh, solve, decideOnRefinement, 10)
  print "Converged with order:", eoc_rec.estimateOrderOfConvergence()[0,1]
  eoc_rec.writeGnuplotFile(",,convergence.data")

  vis.writeMatlabFile("/tmp/visualize.m", new_mesh.dofManager(), new_mesh.elements(), solution_vector)
  #vis.writeGnuplotFile("+result.dat", dof_manager, elements, solution)
  #vis.writeVtkFile("+result.vtk", dof_manager, elements, solution)




adaptiveDemo()

def doProfile():
  profile.run("poissonDemo()", ",,profile")
  p = pstats.Stats(",,profile")
  p.sort_stats("time").print_stats()


