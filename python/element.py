import math

import Numeric as num
import LinearAlgebra as la

import integration
import expression




class tFiniteElementError(Exception):
  def __init__( self, value ):
    self.value = value
  def __str__( self ):
    return repr( self.value )




# node ------------------------------------------------------------------------
class tNode:
  def __init__( self, coordinates ):
    self.Coordinates = coordinates

  def coordinates( self ):
    return self.Coordinates




# finite element abstract interface -------------------------------------------
class tFiniteElement:
  """This class is the abstract interface for an element. An element only
  knows how to add various integral contributions to a matrix
  given indirectly through a tMatrixBuilder."""
  def __init__( self, nodes, dof_manager ):
    self.Nodes = nodes
    self.NodeNumbers = map( dof_manager.getDegreeOfFreedomNumber, self.Nodes )

  def nodes( self ):
    return self.Nodes

  def nodeNumbers( self ):
    return self.NodeNumbers

  def addVolumeIntegralOverDifferentiatedFormFunctions( self, builder, which_derivative = "both" ):
    """This functions adds to the matrix built by `builder' the term 

    \int_{Element} d/dx \phi_i(x,y) d/dx \phi_j(x,y) d(x,y) (for which_derivative == "x")
    \int_{Element} d/dy \phi_i(x,y) d/dy \phi_j(x,y) d(x,y) (for which_derivative == "y")

    where \phi_i and \phi_j run through all the form functions present in
    the element. The correct entries in the matrix are found through the
    DOF manager lookup facility.
    """
    pass

  def addVolumeIntegralOverFormFunctions( self, builder ):
    """This functions adds to the matrix built by `builder' the term 

    \int_{Element} \phi_i(x,y) \phi_j(x,y) d(x,y)

    where \phi_i and \phi_j run through all the form functions present in
    the element. The correct entries in the matrix are found through the
    DOF manager lookup facility.
    """
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
    DOF manager lookup facility.
    """
    pass

  def getSolutionFunction( self, solution_vector ):
    """Once the linear system has been solved, you can use this
    function to obtain the actual solution function which, in
    general, would be the appropriate linear combination of its
    form functions.
    """
    pass




# implementations -------------------------------------------------------------
class tTwoDimensionalLinearTriangularFiniteElement( tFiniteElement ):
  # form function compilation -------------------------------------------------
  FormFunctionExpressions = [ 
    ("-",1,("+",("variable","x"),("variable","y"))),
    ("variable","x"),
    ("variable","y"),
    ]

  FormFunctionCount = len( FormFunctionExpressions )

  def ffcompile( expr ):
    return expression.compile( expr, { "x": "point[0]", "y": "point[1]" }, [ "point" ] )

  FormFunctions = \
    [ ffcompile( expr ) for expr in FormFunctionExpressions ]
  FormFunctionsDx = \
    [ ffcompile( expression.simplify( expression.differentiate( expr, "x" ) ) )
      for expr in FormFunctionExpressions ]
  FormFunctionsDy = \
    [ ffcompile( expression.simplify( expression.differentiate( expr, "y" ) ) )
      for expr in FormFunctionExpressions ]

  def computeFormFunctionCrossIntegrals( formfunc ):
    result = num.zeros( 
	( len( formfunc), len( formfunc )),
	num.Float )

    for i in range( 0, len( formfunc ) ):
      for j in range( 0, i+1 ):
	fi = formfunc[ i ]
	fj = formfunc[ j ]

	result[i,j] = result[j,i] = \
	  integration.integrateFunctionOnUnitTriangle(
	      lambda point: fi( point ) * fj( point ) )
    return result

  FormFunctionCrossIntegrals = computeFormFunctionCrossIntegrals( 
      FormFunctions )

  # initialization ------------------------------------------------------------
  def __init__( self, nodes, dof_manager ):
    tFiniteElement.__init__( self, nodes, dof_manager )

    self.X = map( lambda node: node.coordinates()[ 0 ], self.Nodes )
    self.Y = map( lambda node: node.coordinates()[ 1 ], self.Nodes )

    x = self.X
    y = self.Y

    self.TransformMatrix = num.array( [
	[ x[1] - x[0], x[2] - x[0] ],
	[ y[1] - y[0], y[2] - y[0] ] ] )
    self.TransformMatrixInverse = la.inverse( self.TransformMatrix )
    self.Origin = nodes[0].coordinates()

    self.Area = 0.5 * math.fabs( la.determinant( self.TransformMatrix ) )
    
  # internal helpers ----------------------------------------------------------
  def transformToUnit( self, point ):
    return num.matrixmultiply( self.TransformMatrixInverse, point - self.Origin )

  # external tools ------------------------------------------------------------
  def barycenter( self ):
    return reduce(
      lambda sum,node: sum+node.coordinates(),
      self.Nodes[1:],
      self.Nodes[0].coordinates() ) * ( 1./len( self.Nodes ) )

  def boundingBox( self ):
    coords = [ node.coordinates() for node in self.Nodes ]
    return num.minimum.reduce( coords ), num.maximum.reduce( coords )

  def isInElement( self, point ):
    # convert point to barycentric
    point_minus_origin = point - self.Nodes[0].coordinates()
    barycentric = num.matrixmultiply( self.TransformMatrixInverse, point_minus_origin )
    # check if point is in
    return num.alltrue( barycentric >= 0 ) and sum( barycentric ) < 1

  # tFiniteElement interface --------------------------------------------------
  def addVolumeIntegralOverDifferentiatedFormFunctions( self, builder, which_derivative = "both" ):
    include_x = 0
    include_y = 0

    a00 = self.TransformMatrix[0,0]
    a01 = self.TransformMatrix[0,1]
    a10 = self.TransformMatrix[1,0]
    a11 = self.TransformMatrix[1,1]

    fdxr = None;fdyr = None;fdxc = None;fdyc = None;

    jacobian = 1/(2*self.Area)

    if which_derivative == "both":
      def functionInIntegral( point ):
	return jacobian * ( \
	  ( a00 * fdxr( point ) + a01 * fdyr( point ) ) * \
	  ( a00 * fdxc( point ) + a01 * fdyc( point ) ) + \
	  ( a10 * fdxr( point ) + a11 * fdyr( point ) ) * \
	  ( a10 * fdxc( point ) + a11 * fdyc( point ) ) 
	  )
    else:
      raise tFiniteElementError, "which_derivative != 'both' not implemented"

    node_count = len( self.Nodes )
    influence_matrix = num.zeros( (node_count, node_count), num.Float )
    for row in range( 0, node_count ):
      for column in range( 0, row + 1 ):
	fdxr = self.FormFunctionsDx[ row ]
	fdxc = self.FormFunctionsDx[ column ]
	fdyr = self.FormFunctionsDy[ row ]
	fdyc = self.FormFunctionsDy[ column ]

	influence_matrix[row,column] = \
	influence_matrix[column,row] = \
	    integration.integrateFunctionOnUnitTriangle( functionInIntegral )

    builder.add( influence_matrix, self.NodeNumbers )

  def addVolumeIntegralOverFormFunctions( self, builder ):
    builder.add( self.Area * self.FormFunctionCrossIntegrals, 
	self.NodeNumbers )

  def addVolumeIntegralOverFormFunction( self, builder, f ):
    # FIXME: this is way past inexact
    # is real numerical integration worth the fuss here?
    influences = num.ones( len( self.Nodes ) ) * self.Area * (1/3.) * f( self.barycenter() )
    builder.add( influences, self.NodeNumbers )

  def getSolutionFunction( self, solution_vector ):
    node_values = num.take( solution_vector, self.NodeNumbers )
    def f( point ):
      result = 0
      for i in range( 0, self.FormFunctionCount ):
	result += self.FormFunctions[ i ]( point ) * node_values[ i ]
      return result
    return f





