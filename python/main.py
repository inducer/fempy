#! /usr/bin/python

import Numeric as num
import LinearAlgebra as la
import spmatrix as sparse
import itsolvers
import math
import pyvtk
import precon




# tools -----------------------------------------------------------------------
def flatten( list ):
  result = []
  for i in list:
    result.extend( i )
  return result




def writeSymmetricMatrixAsCSV( filename, matrix ):
  mat_file = file( filename, "w" )
  w,h = matrix.shape
  for row in range( 0, h ):
    for column in range( 0, row ):
      mat_file.write( "%f," % matrix[ row, column ] )
    for column in range( row, w ):
      mat_file.write( "%f," % matrix[ column, row ] )
    mat_file.write( "\n" )




def norm2( vector ):
  return math.sqrt( num.dot( vector, vector ) )




# matrix building -------------------------------------------------------------
class tDOFManager:
  def __init__( self ):
    self.DegreesOfFreedom = { }

  def obtainDegreeOfFreedom( self, identifier ):
    if identifier in self.DegreesOfFreedom:
      return self.DegreesOfFreedom[ identifier ]
    else:
      new_dof_id = len( self.DegreesOfFreedom )
      self.DegreesOfFreedom[ identifier ] = new_dof_id
      return new_dof_id

  def countDegreesOfFreedom( self ):
    return len( self.DegreesOfFreedom )




class tMatrixBuilder:
  def add( self, small_matrix, small_matrix_rows, small_matrix_columns = None ):
    pass

  def matrix( self ):
    """The result of this function is, in general, of unspecified type.
    However, the returned object will at least support the following
    subset of the numarray interface:

      m.shape
      subscripting and slicing
    """
    pass

  def column( self, i ):
    """Returns the i'th column of the built matrix as a copied dense array."""
    pass

  def forceIdentityMap( self, dof_number ):
    mat = self.Matrix
    w,h = mat.shape

    for i in range( 0, w ):
      mat[ dof_number, i ] = 0
    for i in range( 0, h ):
      mat[ i, dof_number ] = 0
    mat[ dof_number, dof_number ] = 1.




class tSymmetricSparseMatrixBuilder( tMatrixBuilder ):
  def __init__( self, size ):
    self.Matrix = sparse.ll_mat_sym( size )

  def matrix( self ):
    return self.Matrix

  def forceIdentityMap( self, dof_number ):
    mat = self.Matrix
    w,h = mat.shape

    # FIXME: optimize with slicing syntax
    for i in range( 0, dof_number ):
      mat[ dof_number, i ] = 0
    for i in range( dof_number + 1, h ):
      mat[ i, dof_number ] = 0
    mat[ dof_number, dof_number ] = 1.

  def column( self, i ):
    w,h = self.Matrix.shape
    col = num.zeros( (h,), num.Float )
    for j in range( 0, i ):
      col[ j ] = self.Matrix[ i,j ]
    for j in range( i, h ):
      col[ j ] = self.Matrix[ j,i ]
    return col

  def add( self, small_matrix, small_matrix_rows, small_matrix_columns = None ):
    self.Matrix.update_add_mask_sym( 
	small_matrix, 
	num.array( small_matrix_rows ),
	num.ones( (len( small_matrix_rows ),) ) )




class tGeneralSparseMatrixBuilder( tMatrixBuilder ):
  def __init__( self, height, width ):
    self.Matrix = sparse.ll_mat( height, width )

  def matrix( self ):
    return self.Matrix

  def add( self, small_matrix, small_matrix_rows, small_matrix_columns = None ):
    if small_matrix_columns is None:
      small_matrix_columns = small_matrix_rows
    self.Matrix.update_add_mask( small_matrix, num.array( small_matrix_rows ), num.array( small_matrix_columns ) )




class tDenseVectorBuilder( tMatrixBuilder ):
  def __init__( self, size ):
    self.Matrix = num.zeros( (size,), num.Float )

  def matrix( self ):
    return self.Matrix

  def add( self, small_matrix, small_matrix_rows, small_matrix_columns = None ):
    for i in range( 0, len( small_matrix_rows ) ):
      if small_matrix[ i ] != 0:
	self.Matrix[ small_matrix_rows[ i ] ] += small_matrix[ i ]




class tDenseMatrixBuilder( tMatrixBuilder ):
  def __init__( self, height, width ):
    self.Matrix = num.zeros( (height,width), num.Float )

  def matrix( self ):
    return self.Matrix

  def add( self, small_matrix, small_matrix_rows, small_matrix_columns = None ):
    w,h = self.Matrix.shape
    if small_matrix_columns is None:
      small_matrix_columns = small_matrix_rows

    for i in range( 0, len( small_matrix_rows ) ):
      for j in range( 0, len( small_matrix_columns ) ):
	self.Matrix[ small_matrix_rows[ i ], small_matrix_columns[ j ] ] += small_matrix[ i ][ j ]

  def forceIdentityMap( self, dof_number ):
    mat = self.Matrix
    mat[ :,dof_number ] = 0
    mat[ dof_number ] = 0
    mat[ dof_number, dof_number ] = 1.




# geometry classes ------------------------------------------------------------
class tNode:
  def __init__( self, coordinates ):
    self.Coordinates = coordinates

  def coordinates( self ):
    return self.Coordinates




# finite element classes ------------------------------------------------------
class tFiniteElementError(Exception):
  def __init__(self, value):
    self.value = value
  def __str__(self):
    return repr(self.value)




class tFiniteElement:
  """This class is the abstract interface for an element. An element only
  knows how to add various integral contributions to a matrix
  given indirectly through a tMatrixBuilder."""
  def __init__( self, nodes, dof_manager ):
    self.Nodes = nodes

  def addVolumeIntegralOverDifferentiatedFormFunctions( self, builder, which_derivative = "both" ):
    """This functions adds to the matrix built by `builder' the term 

    \int_{Element} d/dx \phi_i(x,y) d/dx \phi_j(x,y) d(x,y) (for which_derivative == "x")
    \int_{Element} d/dy \phi_i(x,y) d/dy \phi_j(x,y) d(x,y) (for which_derivative == "y")

    where \phi_i and \phi_j run through all the form functions present in
    the element. The correct entries in the matrix are found through the
    builder's lookup facility."""
    pass

  def addVolumeIntegralOverFormFunctions( self, builder ):
    """This functions adds to the matrix built by `builder' the term 

    \int_{Element} \phi_i(x,y) \phi_j(x,y) d(x,y)

    where \phi_i and \phi_j run through all the form functions present in
    the element. The correct entries in the matrix are found through the
    builder's lookup facility."""
    pass

  def addVolumeIntegralOverFormFunction( self, builder, f ):
    """This functions adds to the matrix built by `builder' the term 

    \int_{Element} f \phi_i(x,y) d(x,y)

    where \phi_i runs through all the form functions present in
    the element and f is a function that accepts a unspecified
    (i.e. array, numarray or list) of the correct length for
    the point of evaluation and returns a single floating
    point value.

    The correct entries in the matrix are found through the
    builder's lookup facility."""
    pass




class tTwoDimensionalTriangularFiniteElement( tFiniteElement ):
  def __init__( self, nodes, dof_manager ):
    tFiniteElement.__init__( self, nodes, dof_manager )
    self.X = map( lambda node: node.coordinates()[ 0 ], self.Nodes )
    self.Y = map( lambda node: node.coordinates()[ 1 ], self.Nodes )

    x = self.X
    y = self.Y

    self.Area = ( 
    (x[1]-x[0])*(y[2]-y[0])-
    (x[2]-x[0])*(y[1]-y[0]) 
    ) * 0.5

  def area( self ):
    return self.Area

  def barycenter( self ):
    return reduce(
      lambda sum,node: sum+node.coordinates(),
      self.Nodes[1:],
      self.Nodes[0].coordinates() ) * ( 1./len( self.Nodes ) )




class tTwoDimensionalLinearTriangularFiniteElement( tTwoDimensionalTriangularFiniteElement ):
  # cf. \cite{50linesofmatlab} 
  # also see related axiom calculation 50_lines_of_matlab_matrix_m.tm

  def __init__( self, nodes, dof_manager ):
    tTwoDimensionalTriangularFiniteElement.__init__( self, nodes, dof_manager )
    x = self.X
    y = self.Y
    mymat2 = num.array( [x,y] )
    mymat = num.concatenate(
       num.ones( (1,3) ), mymat2 )
    self.InvCoordinateMatrix = la.inverse( num.concatenate(
       ( num.ones( (1,3) ), num.array( [x,y] ) ) ) )

    self.NodeIds = map( dof_manager.obtainDegreeOfFreedom, self.Nodes )

  def addVolumeIntegralOverDifferentiatedFormFunctions( self, builder, which_derivative = "both" ):
    include_x = 0
    include_y = 0

    if which_derivative == "x":
      include_x = 1
    elif  which_derivative == "y":
      include_y = 1
    elif which_derivative == "both":
      include_x = 1
      include_y = 1
    else:
      raise tFiniteElementError, "addVolumeIntegral: which_derivative invalid"

    g = num.matrixmultiply( self.InvCoordinateMatrix, num.array(
      [ [ 0, 0 ], [ include_y, 0 ], [ 0, include_x ] ] ) )
    influence_matrix = self.Area * \
      num.matrixmultiply( g, num.transpose( g ) )
    builder.add( influence_matrix, self.NodeIds )

  def addVolumeIntegralOverFormFunctions( self, builder ):
    raise tFiniteElementError, "NYI: addVolumeIntegralOverFormFunctions" 

  def addVolumeIntegralOverFormFunction( self, builder, f ):
    # FIXME: this is way past inexact
    # is real numerical integration worth the fuss here?
    influences = num.ones( len( self.Nodes ) ) * self.Area * (1/3.) * f( self.barycenter() )
    builder.add( influences, self.NodeIds )





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

  print "grad..."
  for el in elements:
    el.addVolumeIntegralOverDifferentiatedFormFunctions( s_builder )
    el.addVolumeIntegralOverFormFunction( b_builder, f )

  b_mat = b_builder.matrix()
  
  print "bcs..."
  for node in dirichlet_nodes:
    boundary_value = u_d( node.coordinates() )
    i = dof_manager.obtainDegreeOfFreedom( node )
    b_mat += s_builder.column( i ) * boundary_value
    s_builder.forceIdentityMap( i )
    b_mat[ i ] = -boundary_value

  negated_b = b_builder.matrix() * -1
  compiled_s = s_builder.matrix().to_sss()
  x = num.zeros( (dof_count,), num.Float )

  print "solving..."
  info, iter, relres = itsolvers.pcg( compiled_s, negated_b, x, 1e-12, dof_count )
  print "info:", info
  print "iter:", iter
  print "relres:", relres

  residual = num.zeros( x.shape, num.Float )
  compiled_s.matvec( x, residual )

  residual -= b_mat
  print norm2( residual )

  return x





# visualization ---------------------------------------------------------------
def writeGnuplotFile( name, nodes, dof_manager, solution ):
  gnuplot_file = file( name, "w" )
  for node in flatten( nodes ):
    gnuplot_file.write( "%f %f %f\n" % (
	node.coordinates()[0],
	node.coordinates()[1],
	solution[ dof_manager.obtainDegreeOfFreedom( node ) ] ) )

def writeVtkFile( name, nodes, values ):
  pass


# driver ----------------------------------------------------------------------
def poissonDemo():
  width = 1.
  height = 1.
  
  nx = 20
  ny = 20

  center = num.array( [ width/2, height/2 ] )

  def f( coord ):
    if norm2( coord - center ) < 0.2:
      return 20
    else:
      return 0

  def u_d( coord ):
    return norm2(coord-center)

  dof_manager = tDOFManager()
  print "geometry..."
  nodes, elements = buildRectangularGeometry( dof_manager, width / nx, height / ny, nx, ny )

  # make the corner nodes dirichlet nodes
  dirichlet_nodes = []
  dirichlet_nodes.extend( nodes[0] )
  dirichlet_nodes.extend( nodes[-1] )
  dirichlet_nodes.extend( map( lambda node_line: node_line[0], nodes[1:-1] ) )
  dirichlet_nodes.extend( map( lambda node_line: node_line[-1], nodes[1:-1] ) )
  
  solution = solvePoisson( dof_manager, elements, dirichlet_nodes, f, u_d )

  writeGnuplotFile( "result.dat", nodes, dof_manager, solution )
  



poissonDemo()
