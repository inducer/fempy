import math

import Numeric as num
import LinearAlgebra as la

import integration
import expression
import form_function
import visualization




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

  def visualizationData( self, solution_vector ):
    """This function returns a visualization.tVisualizationData structure
    for this element, taking into account the given solution vector.
    """
    pass




# implementations -------------------------------------------------------------
class tTwoDimensionalTriangularFiniteElement( tFiniteElement ):
  # form function compilation -------------------------------------------------
  # defined by subclass:
  FormFunctionExpressions = None

  def compileFormFunctions( cls ):
    n = cls.FormFunctionCount = len( cls.FormFunctionExpressions )

    def ffcompile( expr ):
      return expression.compile( expr, { "0": "point[0]", "1": "point[1]" }, [ "point" ] )

    cls.FormFunctions = \
      [ ffcompile( expr ) for expr in cls.FormFunctionExpressions ]
    cls.FormFunctionsDx = \
      [ ffcompile( expression.simplify( expression.differentiate( expr, "0" ) ) )
        for expr in cls.FormFunctionExpressions ]
    cls.FormFunctionsDy = \
      [ ffcompile( expression.simplify( expression.differentiate( expr, "1" ) ) )
        for expr in cls.FormFunctionExpressions ]

    crossint = num.zeros( (n,n), num.Float )

    for i in range( n ):
      for j in range( 0, i+1 ):
        fi = cls.FormFunctions[ i ]
        fj = cls.FormFunctions[ j ]

        crossint[i,j] = crossint[j,i] = \
          integration.integrateFunctionOnUnitTriangle(
              lambda point: fi( point ) * fj( point ) )

    cls.FormFunctionCrossIntegrals = crossint

  compileFormFunctions = classmethod( compileFormFunctions )

  # initialization ------------------------------------------------------------
  def __init__( self, nodes, dof_manager ):
    tFiniteElement.__init__( self, nodes, dof_manager )
    assert len( nodes ) == len( self.FormFunctions )

    self.X = [ node.coordinates()[0] for node in self.Nodes ]
    self.Y = [ node.coordinates()[1] for node in self.Nodes ]

    x = self.X
    y = self.Y

    self.TransformMatrix = num.array( [
	[ x[1] - x[0], x[2] - x[0] ],
	[ y[1] - y[0], y[2] - y[0] ] ] )
    self.TransformMatrixInverse = la.inverse( self.TransformMatrix )
    self.Origin = nodes[0].coordinates()

    self.Area = 0.5 * math.fabs( la.determinant( self.TransformMatrix ) )
    
  # internal helpers ----------------------------------------------------------
  def transformToReal( self, point ):
    return num.matrixmultiply( self.TransformMatrix, point ) + self.Origin
  def transformToUnit( self, point ):
    return num.matrixmultiply( self.TransformMatrixInverse, point - self.Origin )

  # external tools ------------------------------------------------------------
  def barycenter( self ):
    return sum( [ nd.coordinates() for nd in self.Nodes[0:3] ] ) / 3.

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
    n = self.FormFunctionCount
    influences = num.zeros( (n,), num.Float )

    jacobian = self.Area * 2

    for i in range( n ):
      ff = self.FormFunctions[i]
      influences[i] = jacobian * integration.integrateFunctionOnUnitTriangle( 
        lambda point: f( self.transformToReal( point ) ) * ff(  point ) )

    builder.add( influences, self.NodeNumbers )

  def getSolutionFunction( self, solution_vector ):
    node_values = num.take( solution_vector, self.NodeNumbers )
    def f( point ):
      result = 0
      for i in range( 0, self.FormFunctionCount ):
	result += self.FormFunctions[ i ]( point ) * node_values[ i ]
      return result
    return f





class tTwoDimensionalLinearTriangularFiniteElement( tTwoDimensionalTriangularFiniteElement ):
  FormFunctionExpressions = form_function.makeFormFunctions( 1, [ [0,0], [1,0], [0,1] ] )
tTwoDimensionalLinearTriangularFiniteElement.compileFormFunctions()





class tTwoDimensionalQuadraticTriangularFiniteElement( tTwoDimensionalTriangularFiniteElement ):
  FormFunctionExpressions = form_function.makeFormFunctions( 2, 
    [ [0,0], [1,0], [0,1], [0.5,0], [0.5,0.5], [0,0.5] ] )

  def visualizationData( self, solution_vector ):
    return visualization.tVisualizationData( 
	self.Nodes, 
	[ solution_vector[ n ] for n in self.NodeNumbers ],
	[ (0,3,5), (3,1,4), (3,4,5), (5,4,2) ] )

tTwoDimensionalQuadraticTriangularFiniteElement.compileFormFunctions()

