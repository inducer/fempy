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




def solveAdaptively(mesh, solve, iteration_limit = None):
  iteration_count = 0
  solution = None
  while True:
    print "--------------------------------------------------------------------"
    print "ITERATION %d: %d elements, %d nodes" % (iteration_count, len(mesh.elements()), len(mesh.nodes()))
    print "--------------------------------------------------------------------"
    element_needs_refining, solution = solve(mesh, solution)

    iteration_count += 1
    if iteration_limit is not None and iteration_count >= iteration_limit:
      return mesh, solution

    job = tJob("refining...")
    mesh_change = mesh.getRefinement(element_needs_refining)

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



  
def getParallelogram(edge_length = 1, x_skew = 0, y_skew = 0):
  """Returns a list of points (as tuples) that represent
  the parallelogram given by the parameters.
  The order of the points is as follows:

    10
    23
  """
  a = edge_length / 2.
  xs = x_skew / 2.
  ys = y_skew / 2.
  return [(a-xs,-a+ys), (a+xs,a+ys), (-a+xs,a-ys), (-a-xs,-a-ys)]

def getCircle(radius, segments):
  """The returned circle has beginning and end in the positive x
  direction.
  """
  result = []
  h = (2*math.pi) / segments
  inc = -math.pi / 4.
  for i in range(segments):
    result.append((radius * math.cos(inc), radius * math.sin(inc)))
    inc += h
  return result

def getUnitCellGeometry(radius, segments = 50, inner_factor = 0.3):
  return [tShapeSection(getParallelogram(radius), True),
          tShapeSection(getCircle(radius * inner_factor, segments), False)]

def getExactCircle(radius):
  sqrt2_inv = math.sqrt(2)/2 * radius
  r_squared = radius * radius
  return [
    # going counterclockwise.
    #
    #  T
    # L R ^
    #  B 

    # right
    tShapeGuide(0, [-sqrt2_inv, sqrt2_inv],
                ("**",("-",r_squared,("**",("variable","t"),2)),0.5)),
    
    # top
    tShapeGuide(1, [sqrt2_inv, -sqrt2_inv],
                ("**",("-",r_squared,("**",("variable","t"),2)),0.5)),

    # left
    tShapeGuide(0, [sqrt2_inv, -sqrt2_inv],
                ("-",("**",("-",r_squared,("**",("variable","t"),2)),0.5))),

    # bottom
    tShapeGuide(1, [-sqrt2_inv,sqrt2_inv],
                ("-",("**",("-",r_squared,("**",("variable","t"),2)),0.5))),
    ]

def getExactUnitCellGeometry(edge_length, segments = 50, inner_factor = 0.3):
  return [tShapeSection(getParallelogram(edge_length), True),
          tShapeSection(getExactCircle(edge_length * inner_factor), False)]




def adaptiveDemo(expr, mesh, max_iterations = 10):
  # prepare and compile solution functions ------------------------------------
  rhs = expression.simplify(expression.laplace(expr, ["0", "1"]))
  rhs_c = expression.compileScalarField(rhs)

  sol_c = expression.compileScalarField(expr)
  grad_sol = expression.grad(expr, ["0", "1"])
  grad_sol_c = [expression.compileScalarField(expression.simplify(x)) for x in grad_sol]

  # build geometry ------------------------------------------------------------
  job = tJob("geometry")
  
  mesh.generate()
  job.done()

  eoc_rec = eoc.tEOCRecorder()

  it_number = tools.tReference(0)

  def solve(new_mesh, start_solution_vector = None):
    if start_solution_vector is None:
      start_solution_vector = num.zeros((new_mesh.dofManager().countDegreesOfFreedom(),), num.Float)

    solution_vector = solvePoisson(new_mesh, new_mesh.boundaryNodes(),
      rhs_c, sol_c, start_solution_vector)

    job = tJob("error")
    my_estimator = error_estimator.tAnalyticSolutionL2ErrorEstimator(
      new_mesh, solution_vector, sol_c)

    errors = map(my_estimator, new_mesh.elements())
    worst_error = max(errors)
    max_strategy_selection = filter(lambda err: err > 0.5 * worst_error, errors)
    if len(max_strategy_selection) <= 5:
      print "...backed out to quantile strategy..."
      errors.sort()
      max_error = errors[int(len(errors) * 0.9)]
    else:
      max_error = 0.5 * worst_error

    def refine_decision(el):
      if my_estimator(el) >= max_error:
        return 1
      else:
        return 0

    eoc_rec.addDataPoint(len(new_mesh.elements())**0.5, my_estimator.estimateTotalError())
    job.done()

    it_number.set(it_number.get() + 1)
    return refine_decision, solution_vector

  new_mesh, solution_vector = solveAdaptively(mesh, solve, max_iterations)
  print "Converged with order:", eoc_rec.estimateOrderOfConvergence()[0,1]
  eoc_rec.writeGnuplotFile(",,convergence.data")

  #vis.writeMatlabFile("/tmp/visualize.m", new_mesh.dofManager(), new_mesh.elements(), solution_vector)
  #vis.writeGnuplotFile("+result.dat", new_mesh.dofManager(), new_mesh.elements(), solution_vector)
  vis.writeVtkFile("+result.vtk", new_mesh.dofManager(), new_mesh.elements(), solution_vector)




sol = ("sin", ("*", 5, ("*", ("**",("variable","0"),2), ("**",("variable","1"),2))))
mesh = tTwoDimensionalMesh(getExactUnitCellGeometry(edge_length = 2))
#mesh = tTwoDimensionalSimplicalMesh(getParallelogram(2,0.2), [])
adaptiveDemo(sol, mesh, max_iterations = 6)

def doProfile():
  profile.run("poissonDemo()", ",,profile")
  p = pstats.Stats(",,profile")
  p.sort_stats("time").print_stats()


