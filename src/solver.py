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
  #visualization.writeGnuplotSparsityPattern(",,s.gnuplot", s)

  compiled_s = num.asarray(s, s.typecode(), num.SparseExecuteMatrix)
  s_op = algo.makeMatrixOperator(compiled_s)

  return solveSPDSystem(s_op, negated_b)

