from matrix_builder import *
from stopwatch import *
from tools import *
import pylinear.algorithms as algo




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




def solvePoisson(mesh, f, u_d = lambda x: 0, start_vector = None):
  """Solve the Poisson equation

  laplace u = f
  with u = u_d in the dirichlet_nodes.
  """

  dof_manager = mesh.dofManager()
  dof_count = len(dof_manager)

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

  job = tJob("bcs: dirichlet")
  b_mat = b_builder.matrix()
  
  for node in filter(lambda node: node.ConstraintId == "dirichlet",
                     dof_manager.constrainedNodes()):
    boundary_value = u_d(node.Coordinates)
    i = node.Number
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
                               number_of_eigenvalues = 20, 
                               number_of_arnoldi_vectors = 40,
                               tolerance = 1e-10,
                               max_iterations = 0):

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
  
  return zip(results.RitzValues, results.RitzVectors)




class tLaplacianEigenproblemSolver:
  def __init__(self, mesh, f = None, g = lambda x: 1., typecode = num.Float):
    """Solve the eigenproblem of the Laplace operator:
    
      laplace u + f * u = g * lambda * u.
    
    For the nodes marked with the constraint id "dirichlet",
    a (boundary) condition of u = 0 is enforced.
    """

    self._Mesh = mesh

    dof_manager = self._Mesh.dofManager()
    dof_count = len(dof_manager)

    self._SBuilder = tSymmetricSparseMatrixBuilder(dof_count, typecode)
    self._MBuilder = tSymmetricSparseMatrixBuilder(dof_count, typecode)

    job = tJob("stiffness matrix")
    for el in self._Mesh.elements():
      el.addVolumeIntegralOverDifferentiatedFormFunctions(self._SBuilder)
    job.done()

    job = tJob("mass matrix")
    for el in self._Mesh.elements():
      el.addVolumeIntegralOverFormFunctions(self._MBuilder, g)
    job.done()

    if f is not None:
      job = tJob("rhs")
      for el in self._Mesh.elements():
        el.addVolumeIntegralOverFormFunctions(self._SBuilder, f)
      job.done()

    job = tJob("bcs, dirichlet")
    for node in filter(lambda node: node.ConstraintId == "dirichlet",
                       dof_manager.constrainedNodes()):
      i = node.Number
      self._SBuilder.matrix()[i] = 0
      self._MBuilder.matrix()[i] = 0
      self._MBuilder.matrix()[:,i] = 0
    job.done()

  def addPeriodicBoundaryConditions(self, periodicity_nodes):
    """`periodicity_nodes' is a list of list of tuples (node_i, factor).
    For each element ((n_1,alpha_1), ..., (n_k, alpha_k)) of the 
    top-level list, the boundary condition

      alpha_1*u(n_1) = - alpha_2*u(n_2)... - alpha_k*u(n_k) = 0.

    The boundary condition will be implented in n_1's row.

    This method may be called more than once as long as only
    the alpha_i change and the remaining structure of 
    `periodicity_nodes' remains the same.
    """
    job = tJob("bcs, periodic")
    for condition in periodicity_nodes:
      main_node, main_factor = condition[0]

      # clear out s and m
      self._SBuilder.matrix()[main_node.Number] = 0 
      self._SBuilder.matrix()[main_node.Number, main_node.Number] = 1
      self._MBuilder.matrix()[main_node.Number] = 0
      self._MBuilder.matrix()[:,main_node.Number] = 0

      for node, factor in condition[1:]:
        self._SBuilder.matrix()[main_node.Number,node.Number] = factor / main_factor
    job.done()

  def solve(self, sigma, number_of_eigenvalues = 20, number_of_arnoldi_vectors = 40):
    s = self._SBuilder.matrix()
    m = self._MBuilder.matrix()

    return shiftAndInvertEigenproblem(sigma, s, m,
                                        number_of_eigenvalues = number_of_eigenvalues,
                                        number_of_arnoldi_vectors = number_of_arnoldi_vectors)

  def computeEigenpairResidual(self, value, vector):
    return tools.norm2(
      num.matrixmultiply(self._SBuilder.matrix(), vector) 
      - value * num.matrixmultiply(self._MBuilder.matrix(), vector))

