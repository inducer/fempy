import pylinear.matrices as num
import fempy.eoc as eoc
import fempy.expression as expression
import fempy.stopwatch as stopwatch
import fempy.tools as tools
import fempy.solver as solver
import fempy.visualization as visualization
import fempy.error_estimator as error_estimator





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

    job = stopwatch.tJob("refining...")
    mesh_change = mesh.getRefinement(element_needs_refining)

    if mesh_change is None:
      job.done()
      return mesh,solution

    mesh = mesh_change.meshAfter()
    mesh.generate()
    job.done()

    #job = stopwatch.tJob("adapting...")
    #solution = mesh_change.changeSolution(solution)
    #job.done()

    solution = None



  
def visualize(mesh, vector):
  #visualization.writeMatlabFile("/tmp/visualize.m", mesh.dofManager(), mesh.elements(), vector)
  #visualization.writeGnuplotFile(",,result.dat", mesh.dofManager(), mesh.elements(), vector)
  visualization.writeVtkFile(",,result.vtk", mesh.dofManager(), mesh.elements(), vector)



def adaptiveDemo(expr, mesh, max_iterations = 10):
  # prepare and compile solution functions ------------------------------------
  rhs = expression.simplify(expression.laplace(expr, ["0", "1"]))
  rhs_c = expression.compileScalarField(rhs)

  sol_c = expression.compileScalarField(expr)
  grad_sol = expression.grad(expr, ["0", "1"])
  grad_sol_c = expression.compileVectorField(grad_sol)

  # build geometry ------------------------------------------------------------
  job = stopwatch.tJob("geometry")
  mesh.generate()
  job.done()

  eoc_rec = eoc.tEOCRecorder()

  it_number = tools.tReference(0)

  def solve(new_mesh, start_solution_vector = None):
    if start_solution_vector is None:
      start_solution_vector = num.zeros((new_mesh.dofManager().countDegreesOfFreedom(),), num.Float)

    solution_vector = solver.solvePoisson(new_mesh, new_mesh.boundaryNodes(),
      rhs_c, sol_c, start_solution_vector)

    job = stopwatch.tJob("error")
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

    eoc_rec.addDataPoint(len(new_mesh.elements())**0.5, my_estimator.estimateTotalError()**0.5)
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

  visualize(new_mesh, solution_vector)




def doProfile():
  profile.run("poissonDemo()", ",,profile")
  p = pstats.Stats(",,profile")
  p.sort_stats("time").print_stats()

