import pylinear.matrices as num
import fempy.eoc as eoc
import fempy.expression as expression
import fempy.stopwatch as stopwatch
import fempy.tools as tools
import fempy.solver as solver
import fempy.visualization as visualization
import fempy.element_norm as element_norm





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
    print "ITERATION %d: %d elements, %d nodes" % (iteration_count, 
                                                   len(mesh.elements()), 
                                                   len(mesh.dofManager()))
    print "--------------------------------------------------------------------"
    element_needs_refining, solution = solve(mesh, solution)

    iteration_count += 1
    if iteration_limit is not None and iteration_count >= iteration_limit:
      return mesh, solution

    job = stopwatch.tJob("refining...")
    mesh_change = mesh.getRefinement(element_needs_refining)

    if mesh_change is None:
      job.done()
      return mesh, solution

    mesh = mesh_change.meshAfter()
    job.done()

    #job = stopwatch.tJob("adapting...")
    #solution = mesh_change.changeSolution(solution)
    #job.done()

    solution = None



  
def visualize(mesh_func):
  #visualization.visualize("gnuplot", ",,result.data", mesh_func)
  #visualization.visualize("matlab", "/tmp/visualize.m", mesh_func)
  visualization.visualize("vtk", (",,result.vtk", ",,result_grid.vtk"), mesh_func)



def adaptiveDemo(expr, mesh, max_iterations = 10):
  # prepare and compile solution functions ------------------------------------
  rhs = expression.simplify(expression.laplace(expr, ["0", "1"]))
  rhs_c = expression.compileScalarField(rhs)

  sol_c = expression.compileScalarField(expr)
  grad_sol = expression.grad(expr, ["0", "1"])
  grad_sol_c = expression.compileVectorField(grad_sol)

  # build geometry ------------------------------------------------------------

  eoc_rec = eoc.tEOCRecorder()

  it_number = tools.tReference(0)

  def solve(new_mesh, start_solution_vector = None):
    if start_solution_vector is None:
      start_solution_vector = num.zeros((len(new_mesh.dofManager()),), num.Float)

    constraints = solver.getDirichletConstraints(new_mesh, sol_c)
    solution_vector = solver.solvePoisson(new_mesh, rhs_c, constraints,
                                          start_solution_vector)

    job = stopwatch.tJob("error")
    #my_estimator = element_norm.makeEnergyErrorNormSquared(grad_sol_c, solution_vector)
    my_estimator = element_norm.makeL2ErrorNormSquared(sol_c, solution_vector)

    errors = map(my_estimator, new_mesh.elements())
    worst_error = max(errors)
    max_strategy_selection = [err for err in errors if err > 0.5 * worst_error]
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

    eoc_rec.addDataPoint(len(new_mesh.elements())**0.5,
                         tools.sumOver(my_estimator, new_mesh.elements())**0.5)
    job.done()

    it_number.set(it_number.get() + 1)
    return refine_decision, solution_vector

  new_mesh, mesh_function = solveAdaptively(mesh, solve, max_iterations)

  print "-------------------------------------------------------"
  print "EOC overall:", eoc_rec.estimateOrderOfConvergence()[0,1]
  print "EOC Gliding means:"
  gliding_means = eoc_rec.estimateOrderOfConvergence(3)
  gliding_means_iterations,dummy = gliding_means.shape
  for i in range(gliding_means_iterations):
    print "Iteration %d: %f" % (i, gliding_means[i,1])
  print "-------------------------------------------------------"

  eoc_rec.writeGnuplotFile(",,convergence.data")

  visualize(mesh_function)




def doProfile():
  profile.run("poissonDemo()", ",,profile")
  p = pstats.Stats(",,profile")
  p.sort_stats("time").print_stats()

