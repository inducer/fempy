#! /usr/bin/python

import profile
import pstats

# System imports --------------------------------------------------------------
import math
import cmath
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

def getCircle(radius, use_exact = True):
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
                ("**",("-",r_squared,("**",("variable","t"),2)),0.5),
                use_exact_elements = use_exact),
    
    # top
    tShapeGuide(1, [sqrt2_inv, -sqrt2_inv],
                ("**",("-",r_squared,("**",("variable","t"),2)),0.5),
                use_exact_elements = use_exact),

    # left
    tShapeGuide(0, [sqrt2_inv, -sqrt2_inv],
                ("-",("**",("-",r_squared,("**",("variable","t"),2)),0.5)),
                use_exact_elements = use_exact),

    # bottom
    tShapeGuide(1, [-sqrt2_inv,sqrt2_inv],
                ("-",("**",("-",r_squared,("**",("variable","t"),2)),0.5)),
                use_exact_elements = use_exact),
    ]

def getUnitCellGeometry(edge_length, segments = 50, inner_factor = 0.3, use_exact = True):
  return [tShapeSection(getParallelogram(edge_length), True),
          tShapeSection(getCircle(edge_length * inner_factor, use_exact), False)]

def getAnnulusGeometry(outer_radius, inner_radius, use_exact = True):
  return [tShapeSection(getCircle(outer_radius, use_exact), True),
          tShapeSection(getCircle(inner_radius, use_exact), True)]



def adaptiveDemo(expr, mesh, max_iterations = 10):
  # prepare and compile solution functions ------------------------------------
  rhs = expression.simplify(expression.laplace(expr, ["0", "1"]))
  rhs_c = expression.compileScalarField(rhs)

  sol_c = expression.compileScalarField(expr)
  grad_sol = expression.grad(expr, ["0", "1"])
  grad_sol_c = expression.compileVectorField(grad_sol)

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
    my_estimator = error_estimator.tAnalyticSolutionH1ErrorEstimator(
      new_mesh, solution_vector, sol_c, grad_sol_c)

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

  print "-------------------------------------------------------"
  print "EOC overall:", eoc_rec.estimateOrderOfConvergence()[0,1]
  print "EOC Gliding means:"
  gliding_means = eoc_rec.estimateOrderOfConvergence(3)
  gliding_means_iterations,dummy = gliding_means.shape
  for i in range(gliding_means_iterations):
    print "Iteration %d: %f" % (i, gliding_means[i,1])
  print "-------------------------------------------------------"

  eoc_rec.writeGnuplotFile(",,convergence.data")

  #vis.writeMatlabFile("/tmp/visualize.m", new_mesh.dofManager(), new_mesh.elements(), solution_vector)
  #vis.writeGnuplotFile("+result.dat", new_mesh.dofManager(), new_mesh.elements(), solution_vector)
  vis.writeVtkFile("+result.vtk", new_mesh.dofManager(), new_mesh.elements(), solution_vector)




def unitCellDemo(mesh, epsilon, sigma, k):
  job = tJob("geometry")
  mesh.generate()
  job.done()

  bnodes = mesh.boundaryNodes()

  def distanceToLine(start_point, direction, point):
    # Ansatz: start_point + alpha * direction 
    # <start_point + alpha * direction - point, direction> = 0!
    alpha = - num.innerproduct(start_point - point, direction)/tools.norm2squared(direction)
    foot_point = start_point + alpha * direction
    return tools.norm2(point - foot_point), alpha

  job = tJob("periodicity")
  periodic_nodes = []
  thresh = 1e-5
  while bnodes:
    my_node = bnodes[0]

    nodes_with_distances = tools.decorate(
      lambda node: distanceToLine(my_node.coordinates(), k, node.coordinates())[0],
      bnodes[1:])
    close_nodes, far_away_nodes = tools.partition(lambda (node, dist): dist < thresh, nodes_with_distances)

    if close_nodes:
      bnodes = map(lambda (node, dist): node, far_away_nodes)
      tagged_nodes = map(lambda (node, dist): node, close_nodes)
    else:
      nodes_with_distances.sort(lambda (node1, dist1), (node2, dist2): cmp(dist1, dist2))
      tagged_nodes = [nodes_with_distances[0][0]]
      bnodes = map(lambda (node, dist): node, nodes_with_distances[1:])

    for tagged_node in tagged_nodes:
      dist = num.innerproduct(tagged_node.coordinates()-my_node.coordinates(),k)
      periodic_nodes.append((my_node, tagged_node, cmath.exp(1j * dist)))
  job.done()

  connections_file = file(",,connections.data", "w")
  for a,b,factor in periodic_nodes:
    connections_file.write("%f\t%f\n" % (a.coordinates()[0], a.coordinates()[1]))
    connections_file.write("%f\t%f\n" % (b.coordinates()[0], b.coordinates()[1]))
    connections_file.write("\n")

  raw_input()

  solver.solveLaplaceEigenproblem(sigma, mesh, [], periodic_nodes, g = epsilon,
                                  typecode = num.Complex)




def runPoissonDemo():
  sol = ("sin", ("*", 5, ("*", ("**",("variable","0"),2), ("**",("variable","1"),2))))
  mesh = tTwoDimensionalMesh(getUnitCellGeometry(edge_length = 2))
  #mesh = tTwoDimensionalMesh(getAnnulusGeometry(1.5, 0.5, False), hole_starts = [[0,0]])
  adaptiveDemo(sol, mesh, max_iterations = 7)



def runEigenDemo():
  mesh = tTwoDimensionalMesh(getUnitCellGeometry(edge_length = 1, inner_factor = 0.3))
  def epsilon(x):
    if tools.norm2(x) < 0.3:
      return 11
    else:
      return 0

  sigma = 0.9
  
  unitCellDemo(mesh, epsilon, sigma, num.array([1,1]))




runEigenDemo()




def doProfile():
  profile.run("poissonDemo()", ",,profile")
  p = pstats.Stats(",,profile")
  p.sort_stats("time").print_stats()



