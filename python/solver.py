import itsolvers
from matrixbuilder import *
from stopwatch import *
from tools import *




def solvePoisson( dof_manager, elements, dirichlet_nodes, f, u_d = lambda x: 0 ):
  """Solve the Poisson equation

  laplace u = f
  with u = u_d in the dirichlet_nodes.
  """

  dof_count = dof_manager.countDegreesOfFreedom()

  s_builder = tSymmetricSparseMatrixBuilder( dof_count )
  b_builder = tDenseVectorBuilder( dof_count )

  job = tJob( "matrix" )
  for el in elements:
    el.addVolumeIntegralOverDifferentiatedFormFunctions( s_builder )
    el.addVolumeIntegralOverFormFunction( b_builder, 
      lambda x,formfunc_value: f(x) * formfunc_value )
  job.done()

  job = tJob( "bcs" )
  b_mat = b_builder.matrix()
  
  for node in dirichlet_nodes:
    boundary_value = u_d( node.coordinates() )
    i = dof_manager.getDegreeOfFreedomNumber( node )
    b_mat += s_builder.column( i ) * boundary_value
    s_builder.forceIdentityMap( i )
    b_mat[ i ] = -boundary_value

  negated_b = b_builder.matrix() * -1
  compiled_s = s_builder.matrix().to_sss()
  x = num.zeros( (dof_count,), num.Float )
  job.done()

  job = tJob( "solve" )
  info, iter, relres = itsolvers.pcg( compiled_s, negated_b, x, 1e-12, dof_count )
  job.done()
  print "  info:", info
  print "  iter:", iter
  print "  relative residual: ", relres

  residual = num.zeros( x.shape, num.Float )
  compiled_s.matvec( x, residual )
  residual -= b_mat

  print "  absolute residual: ", norm2( residual )

  return x




def solveHelmholtz( dof_manager, elements, dirichlet_nodes, f, u_d = lambda x: 0 ):
  dof_count = dof_manager.countDegreesOfFreedom()

  s_builder = tSymmetricSparseMatrixBuilder( dof_count )
  m_builder = tSymmetricSparseMatrixBuilder( dof_count )

  print "matrix..."
  for el in elements:
    el.addVolumeIntegralOverDifferentiatedFormFunctions( s_builder )
    el.addVolumeIntegralOverFormFunctions( m_builder )









