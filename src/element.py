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

  # Tools ---------------------------------------------------------------------
  def nodes(self):
    return self.Nodes

  def nodeNumbers(self):
    return self.NodeNumbers

  def transformToReal(self):
    """Transforms a unit vector inside the unit element to its corresponding
    vector in the transformed element.
    """
    raise RuntimeError, "not implemented"

  def transformToUnit(self):
    """The inverse of transformToReal.
    """
    raise RuntimeError, "not implemented"

  def getTransformJacobian(self, point):
    """Returns the Jacobian matrix of transformToReal. `point' is in unit 
    coordinates.
    """
    raise RuntimeError, "not implemented"

  def area(self):
    """Returns the area occupied by the transformed element."""
    raise RuntimeError, "not implemented"

  def isInElement(self, point):
    """Returns whether `point' in transformed coordinates is inside the
    element.
    """
    raise RuntimeError, "not implemented"

  # Integral contributions ----------------------------------------------------
  def addVolumeIntegralOverDifferentiatedFormFunctions(self, builder, which_derivative = "both", factor = 1.):
    """This functions adds to the matrix built by `builder' the term 

    \int_{Element} d/dx \phi_i(x,y) d/dx \phi_j(x,y) d(x,y) (for which_derivative == "x")
    \int_{Element} d/dy \phi_i(x,y) d/dy \phi_j(x,y) d(x,y) (for which_derivative == "y")

    where \phi_i and \phi_j run through all the form functions present in
    the element. The correct entries in the matrix are found through the
    DOF manager lookup facility. A sum of both contributions is added if
    which_derivatives is "both".
    """
    raise RuntimeError, "not implemented"

  def addVolumeIntegralOverFormFunctions(self, builder, f):
    """This functions adds to the matrix built by `builder' the term 

    \int_{Element} f(x,y) \phi_i(x,y) \phi_j(x,y) d(x,y)

    where \phi_i and \phi_j run through all the form functions present in
    the element. The correct entries in the matrix are found through the
    DOF manager lookup facility.

    The argument to f is given in transformed coordinates.
    """
    raise RuntimeError, "not implemented"

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
    raise RuntimeError, "not implemented"

  def getVolumeIntegralOver(self, f, coefficients):
    """This functions returns the value of

    \int_{Element} f((x,y), u(x,y)) d(x,y)

    where u is a linear combination of the form functions given
    by the coefficients sequence and f is a function supplied
    by the user.

    (x,y) supplied are in unit coordinates.
    """
    raise RuntimeError, "not implemented"

  # Form functions ------------------------------------------------------------
  def getFormFunctionCombination(self, coefficients):
    """This function returns a function that represents linear combination 
    of form functions, with the given coefficients. The returned functions is 
    defined on the unit element.
    """
    def f(point):
      result = 0
      for i in range(len(self.FormFunctions)):
	result += self.FormFunctions[i](point) * coefficients[ i ]
      return result
    return f

  def getFormFunctionCombinationGradient(self, coefficients):
    """This function returns the gradient of a linear combination of form 
    functions, with the given coefficients. The returned function is defined 
    on the unit element.
    """
    def f(point):
      result = num.zeros((2,), num.Float)
      for i in range(len(self.FormFunctions)):
        for j in range(2):
          result[j] += self.DifferentiatedFormFunctions[i][j](point) * coefficients[i]
      return result
    return f

  # Visualization -------------------------------------------------------------
  def visualizationData(self, solution_vector):
    """This function returns a visualization.tVisualizationData structure
    for this element, taking into account the given solution vector.
    """
    raise RuntimeError, "not implemented"




# Form function related -------------------------------------------------------
def addFormFunctions(cls, form_func_expr, dimensions):
  cls.FormFunctionExpressions = form_func_expr

  cls.FormFunctions = \
    [ expression.compileScalarField(expr) for expr in cls.FormFunctionExpressions ]
  cls.DifferentiatedFormFunctions = \
    [
    [ expression.compileScalarField(expression.simplify(expression.differentiate(expr, "%i" % dim)))
      for dim in range(dimensions) ]
    for expr in cls.FormFunctionExpressions ]




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
    self.InverseDeterminant = 1./determinant
    
  # Tools ---------------------------------------------------------------------
  def transformToReal(self, point):
    return num.matrixmultiply(self.TransformMatrix, point) + self.Origin

  def transformToUnit(self, point):
    return num.matrixmultiply(self.TransformMatrixInverse, point - self.Origin)

  def getTransformJacobian(self, point):
    return self.TransformMatrix

  def area(self):
    return self.Area

  def boundingBox(self):
    coords = [ node.coordinates() for node in self.Nodes ]
    return reduce(num.minimum, coords), reduce(num.maximum, coords)

  def isInElement(self, point):
    unit_coords = self.transformToUnit(point)
    neg_bound = -1e-6
    pos_bound = 1+1e-6
    return \
        neg_bound < unit_coords[0] < pos_bound and \
        neg_bound < unit_coords[1] < pos_bound and \
        neg_bound < 1-unit_coords[0]-unit_coords[1] < pos_bound

  # Internal tools ------------------------------------------------------------
  def _barycenter(self):
    return sum([ nd.coordinates() for nd in self.Nodes[0:3] ]) / 3.

  # Integral contributions ----------------------------------------------------
  def addVolumeIntegralOverDifferentiatedFormFunctions(self, builder, which_derivative = "both", factor = 1.):
    g00 = self.TransformMatrix[0,0]
    g01 = self.TransformMatrix[0,1]
    g10 = self.TransformMatrix[1,0]
    g11 = self.TransformMatrix[1,1]

    if which_derivative == "both":
      def functionInIntegral(point):
	return (\
	  (g11 * fdxr(point) - g10 * fdyr(point)) * \
	  (g11 * fdxc(point) - g10 * fdyc(point)) + \
	  (-g01 * fdxr(point) + g00 * fdyr(point)) * \
	  (-g01 * fdxc(point) + g00 * fdyc(point)) 
	 )
    else:
      raise tFiniteElementError, "which_derivative != 'both' not implemented yet"

    node_count = len(self.Nodes)
    influence_matrix = num.zeros((node_count, node_count), num.Float)
    for row in range(0, node_count):
      for column in range(0, row + 1):
	fdxr = self.DifferentiatedFormFunctions[row][0]
	fdxc = self.DifferentiatedFormFunctions[column][0]
	fdyr = self.DifferentiatedFormFunctions[row][1]
	fdyc = self.DifferentiatedFormFunctions[column][1]

	influence_matrix[row,column] = \
	influence_matrix[column,row] = \
	    integration.integrateOnUnitTriangle(functionInIntegral)

    builder.addScatteredSymmetric(factor * self.InverseDeterminant * influence_matrix, self.NodeNumbers)

  def addVolumeIntegralOverFormFunctions(self, builder, f):
    jacobian_det = self.Area * 2
    def functionInIntegral(point):
      return f(self.transformToReal(point)) * fr(point) * fc(point)

    node_count = len(self.Nodes)
    influence_matrix = num.zeros((node_count, node_count), num.Float)
    for row in range(0, node_count):
      for column in range(0, row + 1):
	fr = self.FormFunctions[row]
	fc = self.FormFunctions[column]

	influence_matrix[row,column] = \
	influence_matrix[column,row] = \
	    integration.integrateOnUnitTriangle(functionInIntegral)

    builder.addScatteredSymmetric(jacobian_det * influence_matrix, 
	self.NodeNumbers)

  def addVolumeIntegralOverFormFunction(self, builder, f):
    n = len(self.FormFunctions)
    influences = num.zeros((n,), num.Float)

    jacobian_det = self.Area * 2

    for i in range(n):
      ff = self.FormFunctions[i]
      influences[i] = jacobian_det * integration.integrateOnUnitTriangle(
        lambda point: f(self.transformToReal(point) , ff( point)))

    builder.addScattered(influences, self.NodeNumbers)

  def getVolumeIntegralOver(self, f, coefficients):
    jacobian_det = self.Area * 2
    zipped = zip(coefficients, self.FormFunctions)
    def functionInIntegral(point):
      ff_comb = sum([ coeff * ff(point) for coeff,ff in zipped])
      return f(point , ff_comb)

    return jacobian_det * integration.integrateOnUnitTriangle(functionInIntegral)




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
  def __init__(self, nodes, distort_function, distort_jacobian, inverse_distort_function, dof_manager):
    tFiniteElement.__init__(self, nodes, dof_manager)
    assert len(nodes) == len(self.FormFunctions)

    self.transformToReal = distort_function
    self.transformToUnit = inverse_distort_function
    self.getTransformJacobian = distort_jacobian

    # verify validity
    if False:
      def inverseNorm(f, inverse_f, point):
        return tools.norm2(point - inverse_f(f(point)))
      inverse_norms = [ 
        inverseNorm(self.transformToReal, self.transformToUnit, num.array(point)) 
        for point in 
          [ [0.,0.], [0.1,0.], [0.1,0.1], [0,0.1], [0.371,0.126], [1.,0.], [0.,1.] ] ]
      print inverse_norms
      assert max(inverse_norms) < 1e-10

  # Tools ---------------------------------------------------------------------
  def area(self):
    def functionInIntegral(point):
      return math.fabs(la.determinant(self.getTransformJacobian(point)))
    return integration.integrateOnUnitTriangle(functionInIntegral)

  def boundingBox(self):
    # FIXME: this will often be wrong
    coords = [ node.coordinates() for node in self.Nodes ]
    return reduce(num.minimum, coords), reduce(num.maximum, coords)

  def isInElement(self, point):
    unit_coords = self.transformToUnit(point)
    neg_bound = -1e-6
    pos_bound = 1+1e-6
    return \
        neg_bound < unit_coords[0] < pos_bound and \
        neg_bound < unit_coords[1] < pos_bound and \
        neg_bound < 1-unit_coords[0]-unit_coords[1] < pos_bound

  # Integral contributions ----------------------------------------------------
  def addVolumeIntegralOverDifferentiatedFormFunctions(self, builder, which_derivative = "both", factor = 1.):
    if which_derivative == "both":
      def functionInIntegral(point):
        g = self.getTransformJacobian(point)

        # determinant count:
        # +1 for the substitution integral
        # -2 (-1 for each occurrence of g, and see: there are two(!) g's multiplied
        #    together in each term)
        # -----------------------------
        # -1.

	return 1/la.determinant(g) * (\
	  (g[1,1] * fdxr(point) - g[1,0] * fdyr(point)) * \
	  (g[1,1] * fdxc(point) - g[1,0] * fdyc(point)) + \
	  (-g[0,1] * fdxr(point) + g[0,0] * fdyr(point)) * \
	  (-g[0,1] * fdxc(point) + g[0,0] * fdyc(point)) 
          )
    else:
      raise tFiniteElementError, "which_derivative != 'both' not implemented yet"

    node_count = len(self.Nodes)
    influence_matrix = num.zeros((node_count, node_count), num.Float)
    for row in range(0, node_count):
      for column in range(0, row + 1):
	fdxr = self.DifferentiatedFormFunctions[row][0]
	fdxc = self.DifferentiatedFormFunctions[column][0]
	fdyr = self.DifferentiatedFormFunctions[row][1]
	fdyc = self.DifferentiatedFormFunctions[column][1]

	influence_matrix[row,column] = \
	influence_matrix[column,row] = \
	    integration.integrateOnUnitTriangle(functionInIntegral)

    builder.addScatteredSymmetric(factor * influence_matrix, self.NodeNumbers)

  def addVolumeIntegralOverFormFunctions(self, builder, f):
    def functionInIntegral(point):
      g = self.getTransformJacobian(point)
      return f(self.transformToReal(point)) * math.fabs(la.determinant(g)) * \
        fr(point) * fc(point)

    node_count = len(self.Nodes)
    influence_matrix = num.zeros((node_count, node_count), num.Float)
    for row in range(0, node_count):
      for column in range(0, row + 1):
	fr = self.FormFunctions[row]
	fc = self.FormFunctions[column]

	influence_matrix[row,column] = \
	influence_matrix[column,row] = \
	    integration.integrateOnUnitTriangle(functionInIntegral)

    builder.addScatteredSymmetric(influence_matrix, self.NodeNumbers)

  def addVolumeIntegralOverFormFunction(self, builder, f):
    n = len(self.FormFunctions)
    influences = num.zeros((n,), num.Float)

    def functionInIntegral(point):
      g = self.getTransformJacobian(point)
      return math.fabs(la.determinant(g)) * \
        f(self.transformToReal(point), ff(point))
    for i in range(n):
      ff = self.FormFunctions[i]
      influences[i] = integration.integrateOnUnitTriangle(functionInIntegral)

    builder.addScattered(influences, self.NodeNumbers)

  def getVolumeIntegralOver(self, f, coefficients):
    zipped = zip(coefficients, self.FormFunctions)
    def functionInIntegral(point):
      ff_comb = sum([ coeff * ff(point) for coeff,ff in zipped])
      return math.fabs(la.determinant(self.getTransformJacobian(point))) * \
          f(point , ff_comb)

    return integration.integrateOnUnitTriangle(functionInIntegral)




# distorted two-dimensional triangular finite elements ------------------------
class tDistortedTwoDimensionalLinearTriangularFiniteElement(tDistortedTwoDimensionalTriangularFiniteElement):
  def visualizationData(self, solution_vector):
    node_numbers_laid_out = []
    nodes = []
    node_values = []

    segments = 8
    h = 1./segments

    for x_n in range(segments+1):
      x = x_n * h
      line_of_node_numbers = []
      for y_n in range(segments-x_n+1):
        y = y_n * h
        nodes.append(self.transformToReal(num.array([x,y])))
        line_of_node_numbers.append(-len(nodes))

        value = 0
        for i in range(len(self.FormFunctions)):
          value += self.FormFunctions[i](num.array([x,y])) * \
                   solution_vector[self.NodeNumbers[i]]
        node_values.append(value)
      node_numbers_laid_out.append(line_of_node_numbers)

    node_numbers_laid_out[0][0] = 0
    node_numbers_laid_out[-1][0] = 1
    node_numbers_laid_out[0][-1] = 2

    triangles = []
    for x_n in range(segments):
      for y_n in range(segments-x_n):
        triangles.append(
          (node_numbers_laid_out[x_n][y_n],
           node_numbers_laid_out[x_n+1][y_n],
           node_numbers_laid_out[x_n][y_n+1]))
        if y_n != segments-x_n-1:
          triangles.append(
            (node_numbers_laid_out[x_n+1][y_n+1],
             node_numbers_laid_out[x_n][y_n+1],
             node_numbers_laid_out[x_n+1][y_n]))

    return visualization.tVisualizationData(
      self.Nodes, 
      [ solution_vector[ n ] for n in self.NodeNumbers ],
      triangles, nodes, node_values)

addFormFunctions(
  tDistortedTwoDimensionalLinearTriangularFiniteElement,
  form_function.makeFormFunctions(1, [ [0,0], [1,0], [0,1] ]),
  dimensions = 2)





