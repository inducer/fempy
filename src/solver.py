from matrix_builder import *
from stopwatch import *
from tools import *
import pylinear.algorithms as algo
import visualization




def solveSPDSystem(matrix_op, rhs, start_vector = None):
  h,w = matrix_op.shape
  if start_vector is None:
    x = num.zeros((w,), num.Float)
  else:
    x = start_vector[:]

  job = tJob("solve")
  matrix_inv_op = algo.makeCGMatrixOperator(matrix_op, h * 2)
  matrix_inv_op.apply(rhs, x)
  job.done()
  print "  iterations:", matrix_inv_op.last_iteration_count

  residual = x.copy()
  matrix_op.apply(x, residual)
  residual -= rhs

  print "  absolute residual: ", norm2(residual)
  return x




def solvePoisson(mesh, dirichlet_nodes, f, u_d = lambda x: 0, start_vector = None):
  """Solve the Poisson equation

  laplace u = f
  with u = u_d in the dirichlet_nodes.
  """

  dof_manager = mesh.dofManager()
  dof_count = dof_manager.countDegreesOfFreedom()

  s_builder = tSymmetricSparseMatrixBuilder(dof_count, num.Float)
  b_builder = tDenseVectorBuilder(dof_count, num.Float)

  job = tJob("matrix")
  for el in mesh.elements():
    el.addVolumeIntegralOverDifferentiatedFormFunctions(s_builder)
  job.done()

  job = tJob("rhs")
  for el in mesh.elements():
    el.addVolumeIntegralOverFormFunction(b_builder, 
      lambda x,formfunc_value: f(x) * formfunc_value)
  job.done()

  job = tJob("bcs")
  b_mat = b_builder.matrix()
  
  for node in dirichlet_nodes:
    boundary_value = u_d(node.coordinates())
    i = dof_manager.getDegreeOfFreedomNumber(node)
    if boundary_value != 0:
      b_mat += s_builder.column(i) * boundary_value
    s_builder.forceIdentityMap(i)
    b_mat[ i ] = -boundary_value

  negated_b = b_builder.matrix() * -1
  job.done()

  s = s_builder.matrix()
  #visualization.writeGnuplotSparsityPattern(",,s.data", s)

  compiled_s = num.asarray(s, s.typecode(), num.SparseExecuteMatrix)
  s_op = algo.makeMatrixOperator(compiled_s)

  return solveSPDSystem(s_op, negated_b)




def shiftAndInvertEigenproblem(sigma, s, m, 
                               number_of_eigenvalues = 5, 
                               number_of_arnoldi_vectors = 10,
                               tolerance = 1e-10,
                               max_iterations = 0):

  compiled_s = num.asarray(s, s.typecode(), num.SparseExecuteMatrix)
  compiled_m = num.asarray(m, m.typecode(), num.SparseExecuteMatrix)

  m_op = algo.makeMatrixOperator(compiled_m)

  job = tJob("shift matrix")
  shifted_matrix = num.asarray(s - sigma * m, s.typecode(), num.SparseExecuteMatrix)
  shifted_matrix_invop = algo.makeUMFPACKMatrixOperator(shifted_matrix)
  job.done()

  op = algo.composeMatrixOperators(shifted_matrix_invop, m_op)

  job = tJob("arpack rci")
  results = algo.runArpack(op, m_op, algo.SHIFT_AND_INVERT_GENERALIZED,
                           sigma, number_of_eigenvalues, number_of_arnoldi_vectors,
                           algo.LARGEST_MAGNITUDE, tolerance, False, max_iterations)
  job.done()
  
  return results




def solveLaplaceEigenproblem(sigma, mesh, 
                             dirichlet_nodes, periodic_nodes = [], 
                             f = None, g = lambda x: 1.,
                             typecode = num.Float):
  """Solve the Poisson equation

  laplace u + f * u = g * lambda * u 

  with u = 0 in nodes given in the list `dirichlet_nodes'. 
  `periodic_nodes' is a list of tuples (node_a, node_b, factor),
  such that a boundary condition of u(node_a) = factor * u(node_b)
  will be forced.
  """

  dof_manager = mesh.dofManager()
  dof_count = dof_manager.countDegreesOfFreedom()

  s_builder = tSymmetricSparseMatrixBuilder(dof_count, typecode)
  m_builder = tSymmetricSparseMatrixBuilder(dof_count, typecode)

  job = tJob("stiffness matrix")
  for el in mesh.elements():
    el.addVolumeIntegralOverDifferentiatedFormFunctions(s_builder)
  job.done()

  job = tJob("mass matrix")
  for el in mesh.elements():
    el.addVolumeIntegralOverFormFunctions(m_builder, g)
  job.done()

  if f is not None:
    job = tJob("mass matrix")
    for el in mesh.elements():
      el.addVolumeIntegralOverFormFunctions(s_builder, f)
    job.done()

  job = tJob("bcs, dirichlet")
  for node in dirichlet_nodes:
    i = dof_manager.getDegreeOfFreedomNumber(node)
    s_builder.forceIdentityMap(i)
    m_builder.matrix()[i] = 0
    m_builder.matrix()[:,i] = 0
  job.done()

  job = tJob("bcs, periodic")
  for node_a, node_b, factor in periodic_nodes:
    nr_a = dof_manager.getDegreeOfFreedomNumber(node_a)
    nr_b = dof_manager.getDegreeOfFreedomNumber(node_b)

    s_builder.forceIdentityMap(nr_a)
    s_builder.matrix()[nr_a,nr_b] = -factor
    m_builder.matrix()[nr_a] = 0
    m_builder.matrix()[:,nr_a] = 0
  job.done()

  s = s_builder.matrix()
  m = m_builder.matrix()

  return shiftAndInvertEigenproblem(sigma, s, m)
