import math

import pylinear.matrices as num
import pylinear.linear_algebra as la

import integration
import expression
import form_function
import visualization
import tools
import mesh_function

class tFiniteElementError(Exception):
    def __init__(self, value):
        Exception.__init__(self, value)





# tNode -----------------------------------------------------------------------
class tNode(object):
    def __init__(self, tag, coordinates = None, 
                 tracking_id = None, shape_section = None):
        self.Coordinates = coordinates
        self.Tag = tag
        self.TrackingId = tracking_id
        self.ShapeSection = shape_section




# tDOFManager -----------------------------------------------------------------
class tDOFManager(object):
    def __init__(self):
        self._TagToNode = {}
        self.Nodes = []
        self._ConstrainedNodes = []

    def registerNode(self, tag, coordinates = None, tracking_id = None, shape_section = None):
        try:
            return self._TagToNode[tag]
        except KeyError:
            node = tNode(tag, coordinates, tracking_id, shape_section)
            self._TagToNode[tag] = node
            self.Nodes.append(node)
            return node

    def getNodeByTag(self, tag):
        try:
            return self._TagToNode[tag]
        except KeyError:
            raise RuntimeError, "Attempted to get a non-registered degree of freedom"
  
    def __len__(self):
        return len(self.Nodes)

    def __getitem__(self, index):
        return self.Nodes[index]




# tools -----------------------------------------------------------------------
def assignNodeNumbers(node_list, beginning_assignment = {}, next_number = None, increment = 1):
    result = beginning_assignment.copy()

    if next_number is None:
        if len(result) == 0:
            next_number = 0
        else:
            next_number = max(result.values())+1

    for node in node_list:
        if not node in result:
            result[node] = next_number
            next_number += increment
    return result




# finite element abstract interface -------------------------------------------
class tFiniteElement(object):
    """This class is the abstract interface for an element. An element only
    knows how to add various integral contributions to a matrix
    given indirectly through a tMatrixBuilder."""

    def __init__(self, base_nodes, dof_manager, form_function_kit):
        self.FormFunctionKit = form_function_kit
        self.FormFunctions = form_function_kit.formFunctions()
        self.DifferentiatedFormFunctions = form_function_kit.differentiatedFormFunctions()
        self.Nodes = form_function_kit.getNodes(base_nodes, dof_manager, self)
    
    def nodes(self):
        return self.Nodes

    def formFunctionKit(self):
        return self.FormFunctionKit

    # Tools ---------------------------------------------------------------------
    def transformToReal(self, point):
        """Transforms a unit vector inside the unit element to its corresponding
        vector in the transformed element.
        """
        raise NotImplementedError

    def transformToUnit(self, point):
        """The inverse of transformToReal.
        """
        raise NotImplementedError

    def getTransformJacobian(self, point):
        """Returns the Jacobian matrix of transformToReal. `point' is in unit 
        coordinates.
        """
        raise NotImplementedError

    def area(self):
        """Returns the area occupied by the transformed element."""
        raise NotImplementedError

    def isInElement(self, point):
        """Returns whether `point' in transformed coordinates is inside the
        element.
        """
        raise NotImplementedError

    # Integral contributions ----------------------------------------------------
    def getVolumeIntegralsOverDifferentiatedFormFunctions(self, which_derivative, factor = 1.):
        """This functions adds to the matrix built by `builder' the term 
        
        \int_{Element} d/dx \phi_i(x,y) d/dx \phi_j(x,y) d(x,y) (for which_derivative == [0])
        \int_{Element} d/dy \phi_i(x,y) d/dy \phi_j(x,y) d(x,y) (for which_derivative == [1])
        ...

        where \phi_i and \phi_j run through all the form functions present in
        the element. The correct entries in the matrix are found through the
        DOF manager lookup facility. 
        Contributions are added for all numbers which occur in which_derivatives.
        which_derivative defaults to all.
        """
        raise NotImplementedError

    def getVolumeIntegralsOverFormFunctions(self, f):
        """This functions adds to the matrix built by `builder' the term 

        \int_{Element} f(x,y) \phi_i(x,y) \phi_j(x,y) d(x,y)

        where \phi_i and \phi_j run through all the form functions present in
        the element. The correct entries in the matrix are found through the
        DOF manager lookup facility.
        
        The argument to f is given in transformed coordinates.
        """
        raise NotImplementedError

    def getVolumeIntegralsOverFormFunction(self, f, typecode = num.Float):
        """This functions adds to the matrix built by `builder' the term 
        
        \int_{Element} f((x,y), \phi_i(x,y)) d(x,y)
        
        where \phi_i runs through all the form functions present in
        the element and f is a function that accepts a sequence
        (i.e. array, numarray or list) of the correct length for
        the point of evaluation and returns a single floating
        point value.

        The (x,y) are supplied in transformed coordinates.

        The correct entries in the matrix are found through the
        DOF manager lookup facility.
        """
        raise NotImplementedError

    def getVolumeIntegralsOver(self, f, coefficients):
        """This functions returns the value of

        \int_{Element} f((x,y), u(x,y)) d(x,y)

        where u is a linear combination of the form functions given
        by the coefficients sequence and f is a function supplied
        by the user.

        The (x,y) are supplied in unit coordinates.
        """
        raise NotImplementedError

    # Visualization -------------------------------------------------------------
    def getVisualizationData(self, solution_vector, vis_data):
        """This function adds itself to the given visualization.tVisualizationData
        data structure, taking into account the given solution vector.
        """
        raise NotImplementedError




# Form function related -------------------------------------------------------
class tInbetweenPoint:
    def __init__(self, node_a_index, node_b_index, numerator, denominator, tag = None):
        self.NodeAIndex = node_a_index
        self.NodeBIndex = node_b_index
        self.Numerator = numerator
        self.Denominator = denominator
        self.Tag = tag

    def completeTag(self, tag_a, tag_b):
        if self.Tag:
            return (tag_a, tag_b, self.Numerator, self.Denominator, self.Tag)
        else:
            return (tag_a, tag_b, self.Numerator, self.Denominator)


class tFormFunctionKit:
    def __init__(self, order, points, extra_points = [], vis_segments = 1):
        dimension = len(points[0])

        all_points = points[:]
        for inbetween_node in extra_points:
            node_a = points[inbetween_node.NodeAIndex]
            node_b = points[inbetween_node.NodeBIndex]
            all_points.append(node_a + (node_b-node_a) * \
                              float(inbetween_node.Numerator) / \
                              float(inbetween_node.Denominator))

        self._Points = points
        self._ExtraPoints = extra_points
        self._VisualizationSegments = vis_segments
    
        self._FormFunctionExpressions = form_function.makeFormFunctions(order,
                                                                    all_points)

        self._FormFunctions = [ expression.compileScalarField(expr, dimension) 
                                for expr in self._FormFunctionExpressions]
        self._DifferentiatedFormFunctions =  [ [ 
            expression.compileScalarField(expression.simplify(
            expression.differentiate(expr, "%i" % dim)), dimension)
            for dim in range(dimension) ] for expr in self._FormFunctionExpressions]
        
    def formFunctions(self):
        return self._FormFunctions

    def differentiatedFormFunctions(self):
        return self._DifferentiatedFormFunctions

    def getNode(self, base_node_list, dof_manager, element, inbetween_point):
        na = base_node_list[inbetween_point.NodeAIndex]
        nb = base_node_list[inbetween_point.NodeBIndex]
        try:
            return dof_manager.getNodeByTag(inbetween_point.completeTag(nb.Tag, na.Tag))
        except RuntimeError:
            pa = self._Points[inbetween_point.NodeAIndex]
            pb = self._Points[inbetween_point.NodeBIndex]
            coordinates = element.transformToReal(pa + (pb-pa) *
                                                  float(inbetween_point.Numerator) /
                                                  float(inbetween_point.Denominator))

            shape_section = None
            tracking_id = None

            if na.ShapeSection is not None and \
               na.ShapeSection is nb.ShapeSection and \
               na.ShapeSection.containsPoint(coordinates):
                shape_section = na.ShapeSection
                tracking_id = shape_section.TrackingId

            return dof_manager.registerNode(inbetween_point.completeTag(na.Tag, nb.Tag),
                                            coordinates,
                                            tracking_id,
                                            shape_section
                                            )

    def getNodes(self, base_nodes, dof_manager, element):
        result = base_nodes + \
                 [self.getNode(base_nodes, dof_manager, element, inbetween_point) 
                  for inbetween_point in self._ExtraPoints]
        return result

    def recommendNumberOfVisualizationSegments(self):
        return self._VisualizationSegments




_Standard1DNodes = [num.array([0]), num.array([1])]
_Standard2DNodes = [num.array([0,0]), num.array([1,0]), num.array([0,1])]

LinearFormFunctionKit1D = tFormFunctionKit(1, _Standard1DNodes)
LinearFormFunctionKit2DTriangle = tFormFunctionKit(1, _Standard2DNodes)
QuadraticFormFunctionKit2DTriangle = tFormFunctionKit(2, _Standard2DNodes,
                                                      [tInbetweenPoint(0, 1, 1, 2),
                                                       tInbetweenPoint(1, 2, 1, 2),
                                                       tInbetweenPoint(2, 0, 1, 2)],
                                                      vis_segments = 2)
# legacy support
LinearFormFunctionKit = LinearFormFunctionKit2DTriangle
QuadraticFormFunctionKit = QuadraticFormFunctionKit2DTriangle




# implementations -------------------------------------------------------------
class tOneDimensionalFiniteElement(tFiniteElement):
    # initialization ------------------------------------------------------------
    def __init__(self, base_nodes, dof_manager, form_function_kit):
        self.Origin = base_nodes[0].Coordinates
        self.End = base_nodes[1].Coordinates
        self.End = base_nodes[1].Coordinates
        self.RawOrigin = num.array([0], num.Float)
        self.RawEnd = num.array([1], num.Float)
        self.Length = (self.End-self.Origin)[0]
        assert self.Length > 0

        tFiniteElement.__init__(self, base_nodes, dof_manager, form_function_kit)
    
    # Tools ---------------------------------------------------------------------
    def transformToReal(self, point):
        return self.Origin + self.Length * point

    def transformToUnit(self, point):
        return (point-self.Origin) / self.Length

    def getTransformJacobian(self, point):
        return num.array([[self.Length]], num.Float)

    def area(self):
        return self.Length

    def boundingBox(self):
        return [self.Origin, self.Origin + num.array([self.Length], num.Float)]

    def isInElement(self, point):
        return 0 <= self.transformToUnit(point) <= 1

    # Integral contributions ----------------------------------------------------
    def getVolumeIntegralsOverDifferentiatedFormFunctions(self, which_derivative = [0], factor = 1.):
        ff_count = len(self.FormFunctions)
        influence_matrix = num.zeros((ff_count, ff_count), num.Float)
        for row in range(0, ff_count):
            for column in range(0, row + 1):
                drdx = self.DifferentiatedFormFunctions[row][0]
                dcdx = self.DifferentiatedFormFunctions[column][0]
                integral = 1/self.Length * \
                           integration.integrateAlongLine(self.RawOrigin, self.RawEnd,
                                                          lambda p: drdx(p)*dcdx(p))
                influence_matrix[row,column] = influence_matrix[column,row] =  integral
        return influence_matrix
        
    def getVolumeIntegralsOverFormFunctions(self, f):
        ff_count = len(self.FormFunctions)
        influence_matrix = num.zeros((ff_count, ff_count), num.Float)
        for row in range(0, ff_count):
            for column in range(0, row + 1):
                rf = self.FormFunctions[row]
                cf = self.FormFunctions[column]
                integral = self.Length * \
                           integration.integrateAlongLine(self.RawOrigin, self.RawEnd,
                                                          lambda p: 
                                                          f(self.transformToReal(p))
                                                          *rf(p)*cf(p))
                influence_matrix[row,column] = influence_matrix[column,row] = \
                                               integral
        return influence_matrix
                          

    def getVolumeIntegralsOverFormFunction(self, f, typecode = num.Float):
        ff_count = len(self.FormFunctions)
        influence_matrix = num.zeros((ff_count,), typecode)
        for row in range(0, ff_count):
            rf = self.FormFunctions[row]
            integral = self.Length * \
                       integration.integrateAlongLine(self.RawOrigin, self.RawEnd,
                                                      lambda p: f(self.transformToReal(p), rf(p)))
            influence_matrix[row] = integral
        return influence_matrix

    def getVolumeIntegralOver(self, f, coefficients):
        ff_count = len(self.FormFunctions)
        influence_matrix = num.zeros((ff_count,), typecode)
        zipped = zip(coefficients, self.FormFunctions)

        def functionInIntegral(p):
            ff_comb = sum([ coeff * ff(point) for coeff,ff in zipped])
            return f(p, ff_comb)

        for row in range(0, ff_count):
            rf = self.FormFunctions[row]
            integral = self.Length * \
                       integration.integrateAlongLine(self.RawOrigin, self.RawEnd,
                                                      functionInIntegral)
            influence_matrix[row] = integral
        return influence_matrix

    def getVisualizationData(self, solution_mesh_func, vis_data):
        raise NotImplementedError





class tTwoDimensionalTriangularFiniteElement(tFiniteElement):
    # initialization ------------------------------------------------------------
    def __init__(self, base_nodes, dof_manager, form_function_kit):
        x = [node.Coordinates[0] for node in base_nodes]
        y = [node.Coordinates[1] for node in base_nodes]

        self.TransformMatrix = num.array([
            [ x[1] - x[0], x[2] - x[0] ],
            [ y[1] - y[0], y[2] - y[0] ] ])
        self.TransformMatrixInverse = la.inverse(self.TransformMatrix)
        self.Origin = base_nodes[0].Coordinates

        determinant = la.determinant(self.TransformMatrix)
        
        # if not, the triangle has the wrong orientation
        assert determinant > 0

        self.Area = 0.5 * determinant
        self.InverseDeterminant = 1./determinant

        tFiniteElement.__init__(self, base_nodes, dof_manager, form_function_kit)
    
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
        coords = [node.Coordinates for node in self.Nodes]
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
        return sum([ nd.Coordinates for nd in self.Nodes[:3] ]) / 3.

    # Integral contributions ----------------------------------------------------
    def getVolumeIntegralsOverDifferentiatedFormFunctions(self, which_derivative = [0,1], factor = 1.):
        g00 = self.TransformMatrix[0,0]
        g01 = self.TransformMatrix[0,1]
        g10 = self.TransformMatrix[1,0]
        g11 = self.TransformMatrix[1,1]

        def functionInIntegral(point):
            return (g11 * fdxr(point) - g10 * fdyr(point)) * \
                   (g11 * fdxc(point) - g10 * fdyc(point)) + \
                   (-g01 * fdxr(point) + g00 * fdyr(point)) * \
                   (-g01 * fdxc(point) + g00 * fdyc(point)) 
        if which_derivative == [0,1]:
            fii = functionInIntegral
        else:
            raise tFiniteElementError, "which_derivative != [0,1] not implemented yet"

        ff_count = len(self.FormFunctions)
        influence_matrix = num.zeros((ff_count, ff_count), num.Float)
        for row in range(0, ff_count):
            for column in range(0, row + 1):
                fdxr = self.DifferentiatedFormFunctions[row][0]
                fdxc = self.DifferentiatedFormFunctions[column][0]
                fdyr = self.DifferentiatedFormFunctions[row][1]
                fdyc = self.DifferentiatedFormFunctions[column][1]
                
                influence_matrix[row,column] = influence_matrix[column,row] = \
                                               integration.integrateOnUnitTriangle(fii)

        return factor * self.InverseDeterminant * influence_matrix

    def getVolumeIntegralsOverFormFunctions(self, f):
        jacobian_det = self.Area * 2
        def functionInIntegral(point):
            return f(self.transformToReal(point)) * fr(point) * fc(point)

        ff_count = len(self.FormFunctions)
        influence_matrix = num.zeros((ff_count, ff_count), num.Float)
        for row in range(0, ff_count):
            for column in range(0, row + 1):
                fr = self.FormFunctions[row]
                fc = self.FormFunctions[column]

                influence_matrix[row,column] = influence_matrix[column,row] = \
                                               integration.integrateOnUnitTriangle(functionInIntegral)

        return jacobian_det * influence_matrix

    def getVolumeIntegralsOverFormFunction(self, f, typecode = num.Float):
        n = len(self.FormFunctions)
        influences = num.zeros((n,), typecode)
        
        jacobian_det = self.Area * 2

        for i in range(n):
            ff = self.FormFunctions[i]
            influences[i] = jacobian_det * integration.integrateOnUnitTriangle(
                lambda point: f(self.transformToReal(point) , ff(point)))

        return influences

    def getVolumeIntegralOver(self, f, coefficients):
        jacobian_det = self.Area * 2
        zipped = zip(coefficients, self.FormFunctions)
        def functionInIntegral(point):
            ff_comb = sum([ coeff * ff(point) for coeff,ff in zipped])
            return f(point , ff_comb)

        return jacobian_det * integration.integrateOnUnitTriangle(functionInIntegral)

    def getVisualizationData(self, solution_mesh_func, vis_data):
        _visualizeTriangle(vis_data, self, solution_mesh_func, 
                           self.FormFunctionKit.recommendNumberOfVisualizationSegments())





# distorted elements ----------------------------------------------------------
class tDistortedTwoDimensionalTriangularFiniteElement(tFiniteElement):
    # initialization ------------------------------------------------------------
    def __init__(self, base_nodes, dof_manager, form_function_kit, 
                 distort_function, distort_jacobian, inverse_distort_function):
        self.transformToReal = distort_function
        self.transformToUnit = inverse_distort_function
        self.getTransformJacobian = distort_jacobian
        tFiniteElement.__init__(self, base_nodes, dof_manager, form_function_kit)

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
        coords = [node.Coordinates for node in self.Nodes]
        return reduce(num.minimum, coords), reduce(num.maximum, coords)

    def isInElement(self, point):
        try:
            unit_coords = self.transformToUnit(point)
            neg_bound = -1e-6
            pos_bound = 1+1e-6
            return \
                   neg_bound < unit_coords[0] < pos_bound and \
                   neg_bound < unit_coords[1] < pos_bound and \
                   neg_bound < 1-unit_coords[0]-unit_coords[1] < pos_bound
        except ValueError:
            # This can happen for deformed elements, the deformation
            # function may not be defined everywhere. If this happens,
            # we decide that we're not inside the element.
            return False

    # Integral contributions ----------------------------------------------------
    def getVolumeIntegralsOverDifferentiatedFormFunctions(self, which_derivative = [0,1], factor = 1.):
        def functionInIntegral(point):
            g = self.getTransformJacobian(point)
            
            # determinant count:
            # +1 for the substitution integral
            # -2 (-1 for each occurrence of g, and see: there are two(!) g's multiplied
            #    together in each term)
            # -----------------------------
            # -1.

            return 1/la.determinant(g) * \
                   ((g[1,1] * fdxr(point) - g[1,0] * fdyr(point)) * \
                    (g[1,1] * fdxc(point) - g[1,0] * fdyc(point)) + \
                    (-g[0,1] * fdxr(point) + g[0,0] * fdyr(point)) * \
                    (-g[0,1] * fdxc(point) + g[0,0] * fdyc(point)))

        if which_derivative == [0,1]:
            fii = functionInIntegral
        else:
            raise tFiniteElementError, "which_derivative != [0,1] not implemented yet"

        ff_count = len(self.FormFunctions)
        influence_matrix = num.zeros((ff_count, ff_count), num.Float)
        for row in range(0, ff_count):
            for column in range(0, row + 1):
                fdxr = self.DifferentiatedFormFunctions[row][0]
                fdxc = self.DifferentiatedFormFunctions[column][0]
                fdyr = self.DifferentiatedFormFunctions[row][1]
                fdyc = self.DifferentiatedFormFunctions[column][1]

                influence_matrix[row,column] = influence_matrix[column,row] = \
                                               integration.integrateOnUnitTriangle(fii)

        return factor * influence_matrix

    def getVolumeIntegralsOverFormFunctions(self, f):
        def functionInIntegral(point):
            g = self.getTransformJacobian(point)
            return f(self.transformToReal(point)) * math.fabs(la.determinant(g)) * \
                   fr(point) * fc(point)

        ff_count = len(self.FormFunctions)
        influence_matrix = num.zeros((ff_count, ff_count), num.Float)
        for row in range(0, ff_count):
            for column in range(0, row + 1):
                fr = self.FormFunctions[row]
                fc = self.FormFunctions[column]

                influence_matrix[row,column] = influence_matrix[column,row] = \
                                               integration.integrateOnUnitTriangle(functionInIntegral)

        return influence_matrix
    
    def getVolumeIntegralsOverFormFunction(self,  f, typecode = num.Float):
        n = len(self.FormFunctions)
        influences = num.zeros((n,), typecode)
        
        def functionInIntegral(point):
            g = self.getTransformJacobian(point)
            return math.fabs(la.determinant(g)) * \
                   f(self.transformToReal(point), ff(point))
        for i in range(n):
            ff = self.FormFunctions[i]
            influences[i] = integration.integrateOnUnitTriangle(functionInIntegral)

        return influences

    def getVolumeIntegralOver(self, f, coefficients):
        zipped = zip(coefficients, self.FormFunctions)
        def functionInIntegral(point):
            ff_comb = sum([ coeff * ff(point) for coeff,ff in zipped])
            return math.fabs(la.determinant(self.getTransformJacobian(point))) * \
                   f(point , ff_comb)

        return integration.integrateOnUnitTriangle(functionInIntegral)

    def getVisualizationData(self, solution_mesh_func, vis_data):
        segments = max(4, self.FormFunctionKit.recommendNumberOfVisualizationSegments())
        _visualizeTriangle(vis_data, self, solution_mesh_func, segments)




# Tools ----------------------------------------------------------------
def _visualizeTriangle(vis_data, element, solution_mesh_func, 
                       segments = 8):
    form_function_kit = element.formFunctionKit()
    form_funcs = form_function_kit.formFunctions()
    element_nodes = element.nodes()
    h = 1./segments

    # first column
    line_of_node_numbers = []

    nuass = solution_mesh_func.numberAssignment()

    line_of_node_numbers.append(nuass[element_nodes[0]])
    for y_n in range(1, segments):
        pt = num.array([0,y_n*h])
      
        line_of_node_numbers.append(
            vis_data.registerLocalNode(
            element.transformToReal(pt), 
            solution_mesh_func.getValueOnElement(element, pt)))

    line_of_node_numbers.append(nuass[element_nodes[2]])
    node_numbers_laid_out = [line_of_node_numbers]

    # in between
    for x_n in range(1,segments):
        x = x_n * h
        line_of_node_numbers = []
        for y_n in range(segments-x_n+1):
            pt = num.array([x,y_n*h])

            line_of_node_numbers.append(
                vis_data.registerLocalNode(
                element.transformToReal(pt), 
                solution_mesh_func.getValueOnElement(element, pt)))

        node_numbers_laid_out.append(line_of_node_numbers)

    # last column
    node_numbers_laid_out.append([nuass[element_nodes[1]]])

    triangles = []
    for x_n in range(segments):
        for y_n in range(segments-x_n):
            vis_data.addTriangle(
                node_numbers_laid_out[x_n][y_n],
                node_numbers_laid_out[x_n+1][y_n],
                node_numbers_laid_out[x_n][y_n+1])
            if y_n != segments-x_n-1:
                vis_data.addTriangle(
                    node_numbers_laid_out[x_n+1][y_n+1],
                    node_numbers_laid_out[x_n][y_n+1],
                    node_numbers_laid_out[x_n+1][y_n])

    # one ccw traversal of the boundary
    my_poly = []
    for x_n in range(segments+1):
        my_poly.append(node_numbers_laid_out[x_n][0])
    for x_n in range(segments)[-1::-1]:
        my_poly.append(node_numbers_laid_out[x_n][segments-x_n])
    for y_n in range(segments)[-1::-1]:
        my_poly.append(node_numbers_laid_out[0][y_n])

    vis_data.addExtraPolygon(my_poly)
