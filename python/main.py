#! /usr/bin/python

import profile

# System imports --------------------------------------------------------------
import math
import sys

# Numerics imports ------------------------------------------------------------
import Numeric as num
import LinearAlgebra as la

# FEM imports -----------------------------------------------------------------
sys.path.append( "pyangle" )
import pyangle
import visualization
import spatial_btree
from matrixbuilder import *
from element import *
from stopwatch import *
from solver import *
from tools import *




# geometry builder ------------------------------------------------------------
def buildRectangularGeometry( dof_manager, dx, dy, nx, ny, second_order = False ):
  # build nodes
  all_nodes = []
  nodes = []
  for node_y in range( 0, ny + 1 ):
    line_nodes = []
    for node_x in range( 0, nx + 1 ):
      line_nodes.append( tNode( num.array( [ node_x * dx, node_y * dy ] ) ) )
    nodes.append( line_nodes )
    all_nodes += line_nodes

  between_nodes = {}

  def between( node1, node2 ):
    if (node1,node2) in between_nodes:
      return between_nodes[ node1,node2 ]
    else:
      new_node = tNode( (node1.coordinates() + node2.coordinates() ) / 2 )
      all_nodes.append( new_node )
      between_nodes[ node1,node2 ] = \
        between_nodes[ node2,node1 ] = new_node
      return new_node

  # build elements, pay attention to mathematically positive orientation
  elements = []
  for el_y in range( 0, ny ):
    for el_x in range( 0, nx ):
      # d c
      # a b
      a = nodes[el_y][el_x]
      b = nodes[el_y][el_x + 1]
      c = nodes[el_y + 1][el_x + 1]
      d = nodes[el_y + 1][el_x]

      if second_order:
        lower_el = tTwoDimensionalQuadraticTriangularFiniteElement( 
          [ a,b,d, between( a, b ), between( b, d ), between( d, a ) ], 
          dof_manager )
        upper_el = tTwoDimensionalQuadraticTriangularFiniteElement( 
          [ c,d,b, between( c, d ), between( d, b ), between( b, c ) ], 
          dof_manager )
      else:
        if False:
          lower_el = makeDistortedLinearElement( [ a,b,d ], dof_manager )
          upper_el = makeDistortedLinearElement( [ c,d,b ], dof_manager )
        else:
          lower_el = tTwoDimensionalLinearTriangularFiniteElement( [ a,b,d ], dof_manager )
          upper_el = tTwoDimensionalLinearTriangularFiniteElement( [ c,d,b ], dof_manager )

      elements.append( lower_el )
      elements.append( upper_el )
      
  return all_nodes, elements
      
  


def buildShapeGeometry( dof_manager, shape_points, refinement_func, second_order = False ):
  out_p = pyangle.triangulateArea( shape_points, refinement_func = refinement_func )

  pts = out_p.Points
  tris = out_p.Triangles

  nodes = []
  for node in range( pts.size() ):
    nodes.append( tNode( num.array( [ pts.getSub( node, 0 ), pts.getSub( node, 1 ) ] ) ) )

  between_nodes = {}

  def between( node1, node2 ):
    if (node1,node2) in between_nodes:
      return between_nodes[ node1,node2 ]
    else:
      new_node = tNode( (node1.coordinates() + node2.coordinates() ) / 2 )
      nodes.append( new_node )
      between_nodes[ node1,node2 ] = \
        between_nodes[ node2,node1 ] = new_node
      return new_node

  elements = []
  for tri in range( tris.size() ):
    a = nodes[ tris.getSub( tri, 0 ) ]
    b = nodes[ tris.getSub( tri, 1 ) ]
    c = nodes[ tris.getSub( tri, 2 ) ]
    if second_order:
      elements.append( tTwoDimensionalQuadraticTriangularFiniteElement( 
        [ a, b, c, between( a, b ), between( b, c ), between( c,a ) ], 
        dof_manager ) )
    else:
      elements.append( tTwoDimensionalLinearTriangularFiniteElement( [a,b,c], dof_manager ) )

  return nodes,elements




def makeDistortedLinearElement( nodes, dof_manager ):
  node_coords = [ node.coordinates() for node in nodes ]

  mat = num.array( [ n - node_coords[0] for n in node_coords[1:] ] )
  matinv = la.inverse( mat )
  nc0 = node_coords[0]

  vars = [ ("variable","0"), ("variable","1") ]
  inv_vars = [ 
    ("-", ("variable","0"), nc0[0]), 
    ("-", ("variable","1"), nc0[1]) ]

  return tDistortedTwoDimensionalLinearTriangularFiniteElement( 
      nodes, [ 
        ("+", nc0[0], expression.linearCombination( mat[0], vars ) ), 
        ("+", nc0[1], expression.linearCombination( mat[1], vars ) ), 
      ], [
        expression.linearCombination( matinv[0], inv_vars ), 
        expression.linearCombination( matinv[1], inv_vars ), 
      ], dof_manager )




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
  
  nx = 40
  ny = 40

  center = num.array( [ width/2, height/2 ] )

  def f( x ):
    if norm2( x - center ) < 0.3:
      return -20
    else:
      return 20

  def u_d( x ):
    if 0.1 < x[0] < 0.9 and x[1] < 0.5:
      return 1
    else:
     return 0
    
    #return math.sin( 5 * ( x[0] + x[1] ) ) + 1
    #return norm2( x - center )

  dof_manager = tDOFManager()

  job = tJob( "geometry" )
  #nodes, elements = buildRectangularGeometry( dof_manager, width / nx, height / ny, nx, ny, False )

  def needsRefinement( vert_origin, vert_destination, vert_apex, area ):
    return area > 0.001
  shape = [ (0,0), (1,0), (1,1), (0,1) ]
  nodes, elements = buildShapeGeometry( dof_manager, shape, needsRefinement, False )
  job.done()

  job = tJob( "btree" )
  finder = spatial_btree.buildElementFinder( elements )
  job.done()

  # make the edge nodes dirichlet nodes
  def isEdgeNode( node ):
    x = node.coordinates()
    return x[0] in [0,1] or x[1] in [0,1]

  dirichlet_nodes = filter( isEdgeNode, nodes )
  
  solution = solvePoisson( dof_manager, elements, dirichlet_nodes, f, u_d )

  s_f1 = makeSolutionFunction( elements, solution, finder )

  #visualization.writeMatlabFile( "/tmp/visualize.m", dof_manager, elements, solution )
  #visualization.writeGnuplotFile( "+result.dat", dof_manager, elements, solution )
  visualization.writeVtkFile( "+result.vtk", dof_manager, elements, solution )
  



def poissonTest():
  width = 1.
  height = 1.
  
  nx = 40
  ny = 40

  a = 5

  def f( point ):
    x = point[0]
    y = point[1]
    return \
      -4 * a **2 * ( x**2 * y**4 + x**4 * y**2 ) * math.sin( a * x**2 * y**2 ) + \
      2 * a * ( y**2 + x**2 ) * math.cos( a * x**2 * y**2 )

  def solution( point ):
    x = point[0]
    y = point[1]
    return math.sin( a * x**2 * y** 2 )
    
  dof_manager = tDOFManager()
  nodes, elements = buildRectangularGeometry( dof_manager, width / nx, height / ny, nx, ny, True )

  # make the edge nodes dirichlet nodes
  def isEdgeNode( node ):
    x = node.coordinates()
    return x[0] in [0,1] or x[1] in [0,1]

  dirichlet_nodes = filter( isEdgeNode, nodes )
  solution_vector = solvePoisson( dof_manager, elements, dirichlet_nodes, f, solution )

  visualization.writeVtkFile( "+result.vtk", dof_manager, elements, solution_vector )

  def errorFunctionL2( point, solution_func_value ):
    result = ( solution( point ) - solution_func_value ) ** 2
    assert result > 0 
    return result

  error_integral = 0
  for el in elements:
    node_values = num.take( solution_vector, el.nodeNumbers() )
    error_integral += el.getVolumeIntegralOver( errorFunctionL2, node_values )

  print "L_2 error estimate:", error_integral
    


#poissonTest()
poissonDemo()
