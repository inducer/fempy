import pytools
import pytools.stopwatch as stopwatch
import pylinear.array as num

import pymbolic

import fempy.eoc as eoc
import fempy.solver as solver
import fempy.visualization as visualization
import fempy.element_norm as element_norm





def solve_adaptively(mesh, solve, iteration_limit = None):
  iteration_count = 0
  solution = None
  while True:
    print "--------------------------------------------------------------------"
    print "ITERATION %d: %d elements, %d nodes" % (iteration_count, 
                                                   len(mesh.elements()), 
                                                   len(mesh.dof_manager()))
    print "--------------------------------------------------------------------"
    element_needs_refining, solution = solve(mesh, solution)

    iteration_count += 1
    if iteration_limit is not None and iteration_count >= iteration_limit:
      return mesh, solution

    job = stopwatch.Job("refining...")
    mesh_change = mesh.get_refinement(element_needs_refining)

    if mesh_change is None:
      job.done()
      return mesh, solution

    mesh = mesh_change.mesh_after()
    job.done()

    #job = stopwatch.Job("adapting...")
    #solution = mesh_change.changeSolution(solution)
    #job.done()

    solution = None



  
def visualize(mesh_func):
  #visualization.visualize("gnuplot", ",,result.data", mesh_func)
  #visualization.visualize("matlab", "/tmp/visualize.m", mesh_func)
  visualization.visualize("vtk", (",,result.vtk", ",,result_grid.vtk"), mesh_func)



def adaptive_demo(expr, mesh, max_iterations = 10):
  variables = [pymbolic.parse("x[0]"), pymbolic.parse("x[1]")]
  # prepare and compile solution functions ------------------------------------
  rhs = pymbolic.simplify(pymbolic.laplace(expr, variables))
  rhs_c = pymbolic.compile(rhs, variables=[pymbolic.var("x")])

  sol_c = pymbolic.compile(expr, variables=[pymbolic.var("x")])
  grad_sol = pymbolic.grad(expr, variables)
  grad_sol_c = pymbolic.VectorFunction(grad_sol, variables=[pymbolic.var("x")])

  # build geometry ------------------------------------------------------------

  energy_eoc_rec = eoc.EOCRecorder()
  l2_eoc_rec = eoc.EOCRecorder()

  it_number = pytools.Reference(0)

  def solve(new_mesh, start_solution_vector = None):
    if start_solution_vector is None:
      start_solution_vector = num.zeros((len(new_mesh.dof_manager()),), num.Float)

    constraints = solver.get_dirichlet_constraints(new_mesh, sol_c)
    solution_vector = solver.solve_poisson(new_mesh, rhs_c, constraints,
                                           start_solution_vector)

    job = stopwatch.Job("error")
    energy_estimator = element_norm.make_energy_error_norm_squared(grad_sol_c, solution_vector)
    l2_estimator = element_norm.make_l2_error_norm_squared(sol_c, solution_vector)

    errors = map(energy_estimator, new_mesh.elements())
    worst_error = max(errors)
    max_strategy_selection = [err for err in errors if err > 0.5 * worst_error]
    if len(max_strategy_selection) <= 5:
      print "...backed out to quantile strategy..."
      errors.sort()
      max_error = errors[int(len(errors) * 0.9)]
    else:
      max_error = 0.5 * worst_error

    def refine_decision(el):
      if energy_estimator(el) >= max_error:
        return 1
      else:
        return 0

    energy_eoc_rec.add_data_point(len(new_mesh.elements())**0.5,
                                  pytools.sum_over(energy_estimator, new_mesh.elements())**0.5)
    l2_eoc_rec.add_data_point(len(new_mesh.elements())**0.5,
                              pytools.sum_over(l2_estimator, new_mesh.elements())**0.5)
    job.done()

    it_number.set(it_number.get() + 1)
    return refine_decision, solution_vector

  new_mesh, mesh_function = solve_adaptively(mesh, solve, max_iterations)

  print "-------------------------------------------------------"

  print "Energy norm EOC overall:", energy_eoc_rec.estimate_order_of_convergence()[0,1]
  print "Energy norm EOC Gliding means:"
  gliding_means = energy_eoc_rec.estimate_order_of_convergence(3)
  gliding_means_iterations,dummy = gliding_means.shape
  for i in range(gliding_means_iterations):
    print "Iteration %d: %f" % (i, gliding_means[i,1])
  print "-------------------------------------------------------"

  energy_eoc_rec.write_gnuplot_file(",,convergence_energy.data")

  print "L^2 EOC overall:", l2_eoc_rec.estimate_order_of_convergence()[0,1]
  print "L^2 EOC Gliding means:"
  gliding_means = l2_eoc_rec.estimate_order_of_convergence(3)
  gliding_means_iterations,dummy = gliding_means.shape
  for i in range(gliding_means_iterations):
    print "Iteration %d: %f" % (i, gliding_means[i,1])
  print "-------------------------------------------------------"

  l2_eoc_rec.write_gnuplot_file(",,convergence_l2.data")

  visualize(mesh_function)
  return mesh_function




def do_profile():
  profile.run("poissonDemo()", ",,profile")
  p = pstats.Stats(",,profile")
  p.sort_stats("time").print_stats()

