#! /usr/bin/python

import profile

# System imports --------------------------------------------------------------
import math
import sys

# Numerics imports ------------------------------------------------------------
import Numeric as num
import LinearAlgebra as la
import spmatrix as sparse
import itsolvers

# FEM imports -----------------------------------------------------------------
sys.path.append( "pyangle" )
import pyangle
import visualization
import spatial_btree
from matrixbuilder import *
from element import *
from stopwatch import *




# tools -----------------------------------------------------------------------
def flatten( list ):
  result = []
  for i in list:
    result.extend( i )
  return result




def writeSymmetricMatrixAsCSV( filename, matrix ):
  mat_file = file( filename, "w" )
  h,w = matrix.shape
  for row in range( 0, h ):
    for column in range( 0, row ):
      mat_file.write( "%f," % matrix[ row, column ] )
    for column in range( row, w ):
      mat_file.write( "%f," % matrix[ column, row ] )
    mat_file.write( "\n" )




def norm2( vector ):
  return math.sqrt( num.dot( vector, vector ) )




def sequence2list( sequence ):
  result = []
  for i in sequence:
    result.append( i )
  return result




def vector2tuple( vector ):
  size, = vector.shape
  if size == 2:
    return vector[ 0 ], vector[ 1 ], 0
  else:
    return vector[ 0 ], vector[ 1 ], vector[ 2 ]




# geometry builder ------------------------------------------------------------
def buildRectangularGeometry( dof_manager, dx, dy, nx, ny ):
  # build nodes
  nodes = [ ]
  for node_y in range( 0, ny + 1 ):
    line_nodes = [ ]
    for node_x in range( 0, nx + 1 ):
      line_nodes.append( tNode( num.array( [ node_x * dx, node_y * dy ] ) ) )
    nodes.append( line_nodes )

  # build elements, pay attention to mathematically positive orientation
  elements = []
  for el_y in range( 0, ny ):
    for el_x in range( 0, nx ):
      lower_el = tTwoDimensionalLinearTriangularFiniteElement(
        [ 
	nodes[ el_y     ][ el_x     ],
	nodes[ el_y     ][ el_x + 1 ],
	nodes[ el_y + 1 ][ el_x     ] ], dof_manager )
      upper_el = tTwoDimensionalLinearTriangularFiniteElement(
        [ 
	nodes[ el_y + 1 ][ el_x + 1 ],
	nodes[ el_y + 1 ][ el_x     ],
	nodes[ el_y     ][ el_x + 1 ] ], dof_manager )

      elements.append( lower_el )
      elements.append( upper_el )
      
  return nodes, elements
      
  


# equation solvers ------------------------------------------------------------
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
    el.addVolumeIntegralOverFormFunction( b_builder, f )
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








# driver ----------------------------------------------------------------------
def makeSolutionFunction( elements, solution, finder = None ):
  if finder is None:
    finder = spatial_btree.buildElementFinder( elements )
  def f( point ):
    el = finder( point )
    if el is None:
      raise RuntimeError, "Solution not defined at %s" % str( point )
    return el.getSolutionFunction( solution )( point )
  return f




def poissonDemo():
  width = 1.
  height = 1.
  
  nx = 20
  ny = 20

  center = num.array( [ width/2, height/2 ] )

  def f( x ):
    if norm2( x - center ) < 0.3:
      return 0
    else:
      return 0

  def u_d( x ):
    if 0.1 < x[0] < 0.9 and x[1] < 0.5:
      return 1
    else:
     return 0
    
    #return math.sin( 5 * ( x[0] + x[1] ) ) + 1
    #return norm2( x - center )

  dof_manager = tDOFManager()

  job = tJob( "geometry" )
  nodes, elements = buildRectangularGeometry( dof_manager, width / nx, height / ny, nx, ny )
  job.done()

  job = tJob( "btree" )
  finder = spatial_btree.buildElementFinder( elements )
  job.done()

  # make the edge nodes dirichlet nodes
  dirichlet_nodes = []
  dirichlet_nodes.extend( nodes[0] )
  dirichlet_nodes.extend( nodes[-1] )
  dirichlet_nodes.extend( map( lambda node_line: node_line[0], nodes[1:-1] ) )
  dirichlet_nodes.extend( map( lambda node_line: node_line[-1], nodes[1:-1] ) )
  
  solution = solvePoisson( dof_manager, elements, dirichlet_nodes, f, u_d )

  s_f1 = makeSolutionFunction( elements, solution, finder )

  visualization.writeMatlabFile( "/tmp/visualize.m", dof_manager, elements, solution )
  #visualization.writeGnuplotFile( "+result.dat", dof_manager, elements, solution )
  #visualisation.writeVtkFile( "+result.vtk", dof_manager, elements, solution )
  



poissonDemo()
