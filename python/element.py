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




# helpers ---------------------------------------------------------------------
def compileExpression( expr ):
  substitutions = {}
  for i in range( 3 ):
    substitutions[ "%d" % i ] = "point[%d]" % i

  return expression.compile( expr, substitutions, [ "point" ] )




def addFormFunctions( cls, form_func_expr, dimensions ):
  cls.FormFunctionExpressions = form_func_expr
  n = cls.FormFunctionCount = len( cls.FormFunctionExpressions )

  cls.FormFunctions = \
    [ compileExpression( expr ) for expr in cls.FormFunctionExpressions ]
  cls.DifferentiatedFormFunctions = \
    [
    [ compileExpression( expression.simplify( expression.differentiate( expr, "%i" % dim ) ) )
      for expr in cls.FormFunctionExpressions ]
    for dim in range( dimensions ) ]

  crossint = num.zeros( (n,n), num.Float )

  for i in range( n ):
    for j in range( 0, i+1 ):
      fi = cls.FormFunctions[ i ]
      fj = cls.FormFunctions[ j ]

      crossint[i,j] = crossint[j,i] = \
        integration.integrateFunctionOnUnitTriangle(
            lambda point: fi( point ) * fj( point ) )

  cls.FormFunctionCrossIntegrals = crossint




# implementations -------------------------------------------------------------
class tTwoDimensionalTriangularFiniteElement( tFiniteElement ):
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

    determinant = la.determinant( self.TransformMatrix )

    # if not, the triangle has the wrong orientation
    assert determinant > 0

    self.Area = 0.5 * determinant
    
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
    g00 = self.TransformMatrixInverse[0,0]
    g01 = self.TransformMatrixInverse[0,1]
    g10 = self.TransformMatrixInverse[1,0]
    g11 = self.TransformMatrixInverse[1,1]

    fdxr = None;fdyr = None;fdxc = None;fdyc = None;

    if which_derivative == "both":
      def functionInIntegral( point ):
	return ( \
	  ( g00 * fdxr( point ) + g10 * fdyr( point ) ) * \
	  ( g00 * fdxc( point ) + g10 * fdyc( point ) ) + \
	  ( g01 * fdxr( point ) + g11 * fdyr( point ) ) * \
	  ( g01 * fdxc( point ) + g11 * fdyc( point ) ) 
	  )
    else:
      raise tFiniteElementError, "which_derivative != 'both' not implemented yet"

    node_count = len( self.Nodes )
    influence_matrix = num.zeros( (node_count, node_count), num.Float )
    for row in range( 0, node_count ):
      for column in range( 0, row + 1 ):
	fdxr = self.DifferentiatedFormFunctions[0][row]
	fdxc = self.DifferentiatedFormFunctions[0][column]
	fdyr = self.DifferentiatedFormFunctions[1][row]
	fdyc = self.DifferentiatedFormFunctions[1][column]

	influence_matrix[row,column] = \
	influence_matrix[column,row] = \
	    integration.integrateFunctionOnUnitTriangle( functionInIntegral )

    jacobian_det = 2*self.Area

    builder.add( jacobian_det * influence_matrix, self.NodeNumbers )

  def addVolumeIntegralOverFormFunctions( self, builder ):
    jacobian_det = self.Area * 2
    builder.add( jacobian_det * self.FormFunctionCrossIntegrals, 
	self.NodeNumbers )

  def addVolumeIntegralOverFormFunction( self, builder, f ):
    n = self.FormFunctionCount
    influences = num.zeros( (n,), num.Float )

    jacobian_det = self.Area * 2

    for i in range( n ):
      ff = self.FormFunctions[i]
      influences[i] = jacobian_det * integration.integrateFunctionOnUnitTriangle( 
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





# concrete two-dimensional triangular elements --------------------------------
class tTwoDimensionalLinearTriangularFiniteElement( tTwoDimensionalTriangularFiniteElement ):
  def visualizationData( self, solution_vector ):
    return visualization.tVisualizationData( 
	self.Nodes, 
	[ solution_vector[ n ] for n in self.NodeNumbers ],
	[ (0,1,2) ] )

addFormFunctions( 
  tTwoDimensionalLinearTriangularFiniteElement,
  form_function.makeFormFunctions( 1, [ [0,0], [1,0], [0,1] ] ),
  dimensions = 2 )





class tTwoDimensionalQuadraticTriangularFiniteElement( tTwoDimensionalTriangularFiniteElement ):
  def visualizationData( self, solution_vector ):
    return visualization.tVisualizationData( 
	self.Nodes, 
	[ solution_vector[ n ] for n in self.NodeNumbers ],
	[ (0,3,5), (3,1,4), (3,4,5), (5,4,2) ] )

addFormFunctions( 
  tTwoDimensionalQuadraticTriangularFiniteElement,
  form_function.makeFormFunctions( 2, 
    [ [0,0], [1,0], [0,1], [0.5,0], [0.5,0.5], [0,0.5] ] ),
  dimensions = 2 )




# distorted elements ----------------------------------------------------------
class tDistortedTwoDimensionalTriangularFiniteElement( tFiniteElement ):
  # initialization ------------------------------------------------------------
  def __init__( self, nodes, distort_expressions, inverse_distort_expressions ):
    tFiniteElement.__init__( self, nodes, dof_manager )
    assert len( nodes ) == len( self.FormFunctions )

    dimensions = 2
    # create derivatives
    self.DistortionDerivatives = []
    for dim in range( dimensions ):
      self.DistortionDerirvatives.append( [] )
      grad = self.DistortionDerivatives[dim]
      f = distort_expressions[dim]

      for ddim in range( dimensions ):
        pass







    # FIXME: verify validity
    
  # internal helpers ----------------------------------------------------------
  def transformToReal( self, point ):
    return num.matrixmultiply( self.TransformMatrix, point ) + self.Origin
  def transformToUnit( self, point ):
    return num.matrixmultiply( self.TransformMatrixInverse, point - self.Origin )

  # external tools ------------------------------------------------------------
  def boundingBox( self ):
    # FIXME: this might be wrong
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

    jacobian_det = 1/(2*self.Area)

    if which_derivative == "both":
      def functionInIntegral( point ):
	return jacobian_det * ( \
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
    jacobian_det = self.Area * 2
    builder.add( jacobian_det * self.FormFunctionCrossIntegrals, 
	self.NodeNumbers )

  def addVolumeIntegralOverFormFunction( self, builder, f ):
    n = self.FormFunctionCount
    influences = num.zeros( (n,), num.Float )

    jacobian_det = self.Area * 2

    for i in range( n ):
      ff = self.FormFunctions[i]
      influences[i] = jacobian_det * integration.integrateFunctionOnUnitTriangle( 
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




# distorted two-dimensional triangular finite elements ------------------------
class tDistortedTwoDimensionalLinearTriangularFiniteElement( tDistortedTwoDimensionalTriangularFiniteElement ):
  def visualizationData( self, solution_vector ):
    return visualization.tVisualizationData( 
	self.Nodes, 
	[ solution_vector[ n ] for n in self.NodeNumbers ],
	[ (0,1,2) ] )

addFormFunctions( 
  tDistortedTwoDimensionalLinearTriangularFiniteElement,
  form_function.makeFormFunctions( 1, [ [0,0], [1,0], [0,1] ] ),
  dimensions = 2 )
