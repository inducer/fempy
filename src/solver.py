from matrixbuilder import *
from stopwatch import *
from tools import *
import pylinear.algorithms as algo




def solvePoisson(dof_manager, elements, dirichlet_nodes, f, u_d = lambda x: 0):
  """Solve the Poisson equation

  laplace u = f
  with u = u_d in the dirichlet_nodes.
  """

  dof_count = dof_manager.countDegreesOfFreedom()

  s_builder = tSymmetricSparseMatrixBuilder(dof_count)
  b_builder = tDenseVectorBuilder(dof_count)

  job = tJob("matrix")
  for el in elements:
    el.addVolumeIntegralOverDifferentiatedFormFunctions(s_builder)
  job.done()
  job = tJob("rhs")
  for el in elements:
    el.addVolumeIntegralOverFormFunction(b_builder, 
      lambda x,formfunc_value: f(x) * formfunc_value)
  job.done()

  job = tJob("bcs")
  b_mat = b_builder.matrix()
  
  nonzero_diri_nodes = []
  for node in dirichlet_nodes:
    boundary_value = u_d(node.coordinates())
    i = dof_manager.getDegreeOfFreedomNumber(node)
    if boundary_value != 0:
      b_mat += s_builder.column(i) * boundary_value
      nonzero_diri_nodes.append((i,boundary_value))
    s_builder.forceIdentityMap(i)
    b_mat[ i ] = -boundary_value

  negated_b = b_builder.matrix() * -1

  s = s_builder.matrix()

  compiled_s = num.asarray(s, s.typecode(), num.SparseExecuteMatrix)
  s_op = algo.makeMatrixOperator(compiled_s)
  s_inv_op = algo.makeCGMatrixOperator(s_op, dof_count)

  x = num.zeros((dof_count,), num.Float)
  job.done()

  job = tJob("solve")
  s_inv_op.apply(negated_b, x)
  job.done()
  print "  iter:", s_inv_op.last_iteration_count

  residual = num.matrixmultiply(compiled_s, x) - b_mat

  print "  absolute residual: ", norm2(residual)

  return x




def solveHelmholtz(dof_manager, elements, dirichlet_nodes, f, u_d = lambda x: 0):
  dof_count = dof_manager.countDegreesOfFreedom()

  s_builder = tSymmetricSparseMatrixBuilder(dof_count)
  m_builder = tSymmetricSparseMatrixBuilder(dof_count)

  print "matrix..."
  for el in elements:
    el.addVolumeIntegralOverDifferentiatedFormFunctions(s_builder)
    el.addVolumeIntegralOverFormFunctions(m_builder)









