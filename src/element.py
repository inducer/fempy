import math

import pylinear.matrices as num
import pylinear.linear_algebra as la

import integration
import expression
import form_function
import visualization
import tools


class tFiniteElementError(Exception):
  def __init__(self, value):
    self.value = value
  def __str__(self):
    return repr(self.value)




# node ------------------------------------------------------------------------
class tNode:
  def __init__(self, coordinates):
    self.Coordinates = coordinates

  def coordinates(self):
    return self.Coordinates




# finite element abstract interface -------------------------------------------
class tFiniteElement:
  """This class is the abstract interface for an element. An element only
  knows how to add various integral contributions to a matrix
  given indirectly through a tMatrixBuilder."""
  def __init__(self, nodes, dof_manager):
    self.Nodes = nodes
    self.NodeNumbers = map(dof_manager.registerDegreeOfFreedomNumber, self.Nodes)

  def nodes(self):
    return self.Nodes

  def nodeNumbers(self):
    return self.NodeNumbers

  def addVolumeIntegralOverDifferentiatedFormFunctions(self, builder, which_derivative = "both"):
    """This functions adds to the matrix built by `builder' the term 

    \int_{Element} d/dx \phi_i(x,y) d/dx \phi_j(x,y) d(x,y) (for which_derivative == "x")
    \int_{Element} d/dy \phi_i(x,y) d/dy \phi_j(x,y) d(x,y) (for which_derivative == "y")

    where \phi_i and \phi_j run through all the form functions present in
    the element. The correct entries in the matrix are found through the
    DOF manager lookup facility.
    """
    pass

  def addVolumeIntegralOverFormFunctions(self, builder):
    """This functions adds to the matrix built by `builder' the term 

    \int_{Element} \phi_i(x,y) \phi_j(x,y) d(x,y)

    where \phi_i and \phi_j run through all the form functions present in
    the element. The correct entries in the matrix are found through the
    DOF manager lookup facility.
    """
    pass

  def addVolumeIntegralOverFormFunction(self, builder, f):
    """This functions adds to the matrix built by `builder' the term 

    \int_{Element} f((x,y), \phi_i(x,y)) d(x,y)

    where \phi_i runs through all the form functions present in
    the element and f is a function that accepts a unspecified
    (i.e. array, numarray or list) of the correct length for
    the point of evaluation and returns a single floating
    point value.

    The correct entries in the matrix are found through the
    DOF manager lookup facility.
    """
    pass

  def getVolumeIntegralOver(self, f, coefficients):
    """This functions returns the value of

    \int_{Element} f((x,y), u(x,y)) d(x,y)

    where u is a linear combination of the form functions given
    by the coefficients sequence and f is a function supplied
    by the user.
    """

  def getSolutionFunction(self, solution_vector):
    """Once the linear system has been solved, you can use this
    function to obtain the actual solution function which, in
    general, would be the appropriate linear combination of its
    form functions.
    """
    pass

  def visualizationData(self, solution_vector):
    """This function returns a visualization.tVisualizationData structure
    for this element, taking into account the given solution vector.
    """
    pass




# helpers ---------------------------------------------------------------------
def addFormFunctions(cls, form_func_expr, dimensions):
  cls.FormFunctionExpressions = form_func_expr
  n = cls.FormFunctionCount = len(cls.FormFunctionExpressions)

  cls.FormFunctions = \
    [ expression.compileScalarField(expr) for expr in cls.FormFunctionExpressions ]
  cls.DifferentiatedFormFunctions = \
    [
    [ expression.compileScalarField(expression.simplify(expression.differentiate(expr, "%i" % dim)))
      for expr in cls.FormFunctionExpressions ]
    for dim in range(dimensions) ]

  crossint = num.zeros((n,n), num.Float)

  for i in range(n):
    for j in range(0, i+1):
      fi = cls.FormFunctions[ i ]
      fj = cls.FormFunctions[ j ]

      crossint[i,j] = crossint[j,i] = \
        integration.integrateFunctionOnUnitTriangle(
            lambda point: fi(point) * fj(point))

  cls.FormFunctionCrossIntegrals = crossint




# implementations -------------------------------------------------------------
class tTwoDimensionalTriangularFiniteElement(tFiniteElement):
  # initialization ------------------------------------------------------------
  def __init__(self, nodes, dof_manager):
    tFiniteElement.__init__(self, nodes, dof_manager)
    assert len(nodes) == len(self.FormFunctions)

    self.X = [ node.coordinates()[0] for node in self.Nodes ]
    self.Y = [ node.coordinates()[1] for node in self.Nodes ]

    x = self.X
    y = self.Y

    self.TransformMatrix = num.array([
	[ x[1] - x[0], x[2] - x[0] ],
	[ y[1] - y[0], y[2] - y[0] ] ])
    self.TransformMatrixInverse = la.inverse(self.TransformMatrix)
    self.Origin = nodes[0].coordinates()

    determinant = la.determinant(self.TransformMatrix)

    # if not, the triangle has the wrong orientation
    assert determinant > 0

    self.Area = 0.5 * determinant
    
  # internal helpers ----------------------------------------------------------
  def transformToReal(self, point):
    return num.matrixmultiply(self.TransformMatrix, point) + self.Origin
  def transformToUnit(self, point):
    return num.matrixmultiply(self.TransformMatrixInverse, point - self.Origin)

  # external tools ------------------------------------------------------------
  def area(self):
    return self.Area

  def barycenter(self):
    return sum([ nd.coordinates() for nd in self.Nodes[0:3] ]) / 3.

  def boundingBox(self):
    coords = [ node.coordinates() for node in self.Nodes ]
    return reduce(num.minimum, coords), reduce(num.maximum, coords)

  def isInElement(self, point):
    # convert point to barycentric
    barycentric = num.matrixmultiply(self.TransformMatrixInverse, 
      point - self.Nodes[0].coordinates())
    # check if point is in
    for i in barycentric:
      if i < 0:
        return False
    return sum(barycentric) < 1

  # tFiniteElement interface --------------------------------------------------
  def addVolumeIntegralOverDifferentiatedFormFunctions(self, builder, which_derivative = "both"):
    g00 = self.TransformMatrixInverse[0,0]
    g01 = self.TransformMatrixInverse[0,1]
    g10 = self.TransformMatrixInverse[1,0]
    g11 = self.TransformMatrixInverse[1,1]

    if which_derivative == "both":
      def functionInIntegral(point):
	return (\
	  (g00 * fdxr(point) + g10 * fdyr(point)) * \
	  (g00 * fdxc(point) + g10 * fdyc(point)) + \
	  (g01 * fdxr(point) + g11 * fdyr(point)) * \
	  (g01 * fdxc(point) + g11 * fdyc(point)) 
	 )
    else:
      raise tFiniteElementError, "which_derivative != 'both' not implemented yet"

    node_count = len(self.Nodes)
    influence_matrix = num.zeros((node_count, node_count), num.Float)
    for row in range(0, node_count):
      for column in range(0, row + 1):
	fdxr = self.DifferentiatedFormFunctions[0][row]
	fdxc = self.DifferentiatedFormFunctions[0][column]
	fdyr = self.DifferentiatedFormFunctions[1][row]
	fdyc = self.DifferentiatedFormFunctions[1][column]

	influence_matrix[row,column] = \
	influence_matrix[column,row] = \
	    integration.integrateFunctionOnUnitTriangle(functionInIntegral)

    jacobian_det = 2*self.Area

    builder.addScatteredSymmetric(jacobian_det * influence_matrix, self.NodeNumbers)

  def addVolumeIntegralOverFormFunctions(self, builder):
    jacobian_det = self.Area * 2
    builder.addScatteredSymmetric(jacobian_det * self.FormFunctionCrossIntegrals, 
	self.NodeNumbers)

  def addVolumeIntegralOverFormFunction(self, builder, f):
    n = self.FormFunctionCount
    influences = num.zeros((n,), num.Float)

    jacobian_det = self.Area * 2

    for i in range(n):
      ff = self.FormFunctions[i]
      influences[i] = jacobian_det * integration.integrateFunctionOnUnitTriangle(
        lambda point: f(self.transformToReal(point) , ff( point)))

    builder.addScattered(influences, self.NodeNumbers)

  def getVolumeIntegralOver(self, f, coefficients):
    jacobian_det = self.Area * 2
    zipped = zip(coefficients, self.FormFunctions)
    def functionInIntegral(point):
      ff_comb = sum([ coeff * ff(point) for coeff,ff in zipped])
      return f(self.transformToReal(point) , ff_comb)

    return jacobian_det * integration.integrateFunctionOnUnitTriangle(functionInIntegral)

  def getSolutionFunction(self, solution_vector):
    node_values = num.take(solution_vector, self.NodeNumbers)
    def f(point):
      unit_point = self.transformToUnit(point)
      result = 0
      for i in range(0, self.FormFunctionCount):
	result += self.FormFunctions[ i ](unit_point) * node_values[ i ]
      return result
    return f





# concrete two-dimensional triangular elements --------------------------------
class tTwoDimensionalLinearTriangularFiniteElement(tTwoDimensionalTriangularFiniteElement):
  def visualizationData(self, solution_vector):
    return visualization.tVisualizationData(
	self.Nodes, 
	[ solution_vector[ n ] for n in self.NodeNumbers ],
	[ (0,1,2) ])

addFormFunctions(
  tTwoDimensionalLinearTriangularFiniteElement,
  form_function.makeFormFunctions(1, [ [0,0], [1,0], [0,1] ]),
  dimensions = 2)





class tTwoDimensionalQuadraticTriangularFiniteElement(tTwoDimensionalTriangularFiniteElement):
  def visualizationData(self, solution_vector):
    return visualization.tVisualizationData(
	self.Nodes, 
	[ solution_vector[ n ] for n in self.NodeNumbers ],
	[ (0,3,5), (3,1,4), (3,4,5), (5,4,2) ])

addFormFunctions(
  tTwoDimensionalQuadraticTriangularFiniteElement,
  form_function.makeFormFunctions(2, 
    [ [0,0], [1,0], [0,1], [0.5,0], [0.5,0.5], [0,0.5] ]),
  dimensions = 2)




# distorted elements ----------------------------------------------------------
class tDistortedTwoDimensionalTriangularFiniteElement(tFiniteElement):
  # initialization ------------------------------------------------------------
  def __init__(self, nodes, distort_expressions, inverse_distort_expressions, dof_manager):
    tFiniteElement.__init__(self, nodes, dof_manager)
    assert len(nodes) == len(self.FormFunctions)

    dimensions = len(distort_expressions)

    self.transformToReal = expression.assembleMatrixFunction([
      expression.compileScalarField(exp) for exp in distort_expressions ])
    self.transformToUnit = expression.assembleMatrixFunction([
      expression.compileScalarField(exp) for exp in inverse_distort_expressions ])

    # create derivatives
    distortion_derivatives = []
    inverse_distortion_derivatives = []
    for dim in range(dimensions):
      grad = []
      inverse_grad = []
      f = distort_expressions[dim]
      inverse_f = inverse_distort_expressions[dim]

      for ddim in range(dimensions):
        grad.append(
          expression.compileScalarField(
            expression.simplify(expression.differentiate(f, "%d" % ddim))))
        inverse_grad.append(
          expression.compileScalarField(
            expression.simplify(expression.differentiate(inverse_f, "%d" % ddim))))
      distortion_derivatives.append(grad)
      inverse_distortion_derivatives.append(inverse_grad)
    self.getDistortionDerivatives = expression.assembleMatrixFunction(
      distortion_derivatives)
    self.getInverseDistortionDerivatives = expression.assembleMatrixFunction(
      inverse_distortion_derivatives)

    # verify validity
    def inverseNorm(f, inverse_f, point):
      return tools.norm2(point - inverse_f(f(point)))
    inverse_norms = [ 
      inverseNorm(self.transformToReal, self.transformToUnit, num.array(point)) 
      for point in 
        [ [0.,0.], [0.1,0.], [0.1,0.1], [0,0.1], [0.371,0.126], [1.,0.], [0.,1.] ] ]
    assert max(inverse_norms) < 1e-11
    
  # external tools ------------------------------------------------------------
  def boundingBox(self):
    # FIXME: this might be wrong
    coords = [ node.coordinates() for node in self.Nodes ]
    return reduce(num.minimum, coords), reduce(num.maximum, coords)

  def isInElement(self, point):
    # FIXME
    raise RuntimeError, "unimplemented"

  # tFiniteElement interface --------------------------------------------------
  def addVolumeIntegralOverDifferentiatedFormFunctions(self, builder, which_derivative = "both"):
    fdxr = None;fdyr = None;fdxc = None;fdyc = None;

    if which_derivative == "both":
      def functionInIntegral(point):
        g = self.getInverseDistortionDerivatives(self.transformToReal(point))
	return 1. / math.fabs(la.determinant(g)) * \
          (\
	  (g[0,0] * fdxr(point) + g[1,0] * fdyr(point)) * \
	  (g[0,0] * fdxc(point) + g[1,0] * fdyc(point)) + \
	  (g[0,1] * fdxr(point) + g[1,1] * fdyr(point)) * \
	  (g[0,1] * fdxc(point) + g[1,1] * fdyc(point)) 
	 )

    else:
      raise tFiniteElementError, "which_derivative != 'both' not implemented yet"

    node_count = len(self.Nodes)
    influence_matrix = num.zeros((node_count, node_count), num.Float)
    for row in range(0, node_count):
      for column in range(0, row + 1):
	fdxr = self.DifferentiatedFormFunctions[0][row]
	fdxc = self.DifferentiatedFormFunctions[0][column]
	fdyr = self.DifferentiatedFormFunctions[1][row]
	fdyc = self.DifferentiatedFormFunctions[1][column]

	influence_matrix[row,column] = \
	influence_matrix[column,row] = \
	    integration.integrateFunctionOnUnitTriangle(functionInIntegral)

    builder.addScatteredSymmetric(influence_matrix, self.NodeNumbers)

  def addVolumeIntegralOverFormFunctions(self, builder):
    def functionInIntegral(point):
      g = self.getDistortionDerivatives(point)
      return math.fabs(la.determinant(g)) * \
        fr(point) * fc(point)

    node_count = len(self.Nodes)
    influence_matrix = num.zeros((node_count, node_count), num.Float)
    for row in range(0, node_count):
      for column in range(0, row + 1):
	fr = self.FormFunctions[row]
	fc = self.FormFunctions[column]

	influence_matrix[row,column] = \
	influence_matrix[column,row] = \
	    integration.integrateFunctionOnUnitTriangle(functionInIntegral)

    builder.addScattered(influence_matrix, self.NodeNumbers)

  def addVolumeIntegralOverFormFunction(self, builder, f):
    n = self.FormFunctionCount
    influences = num.zeros((n,), num.Float)

    def functionInIntegral(point):
      g = self.getDistortionDerivatives(point)
      return math.fabs(la.determinant(g)) * \
        f(self.transformToReal(point), ff(point))
    for i in range(n):
      ff = self.FormFunctions[i]
      influences[i] = integration.integrateFunctionOnUnitTriangle(functionInIntegral)

    builder.addScattered(influences, self.NodeNumbers)

  def getSolutionFunction(self, solution_vector):
    node_values = num.take(solution_vector, self.NodeNumbers)
    def f(point):
      unit_point = self.transformToUnit(point)
      result = 0
      for i in range(0, self.FormFunctionCount):
	result += self.FormFunctions[ i ](unit_point) * node_values[ i ]
      return result
    return f




# distorted two-dimensional triangular finite elements ------------------------
class tDistortedTwoDimensionalLinearTriangularFiniteElement(tDistortedTwoDimensionalTriangularFiniteElement):
  def visualizationData(self, solution_vector):
    return visualization.tVisualizationData(
	self.Nodes, 
	[ solution_vector[ n ] for n in self.NodeNumbers ],
	[ (0,1,2) ])

addFormFunctions(
  tDistortedTwoDimensionalLinearTriangularFiniteElement,
  form_function.makeFormFunctions(1, [ [0,0], [1,0], [0,1] ]),
  dimensions = 2)

