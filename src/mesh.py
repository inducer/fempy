import math

import pytools

import pyangle

import pymbolic
import pymbolic.vector

import pylinear.array as num
import pylinear.computation as comp
import pylinear.toybox as toybox

import fempy.element as element
import fempy.spatial_btree as spatial_btree




# Mesh description data structures ------------------------------------------------
class ShapeGuide:
    """Represents a non-straight-line segment of a shape section.
    """
    def __init__(self, deformation_coordinate, interval, expression, 
            initial_point_count = 3, render_final_point = False, use_exact_elements = True):
        self.DeformationCoordinate = deformation_coordinate
        self.Interval = interval
        self.Expression = expression
        self.InitialPointCount = initial_point_count
        self.RenderFinalPoint = render_final_point
        self.UseExactElements = use_exact_elements

    def end_point(self):
        a,b = self.Interval
        if self.RenderFinalPoint:
            return self.evaluate_to_vector(b)
        else:
            h = float(b-a)/self.InitialPointCount
            return self.evaluate_to_vector(a+h*(self.InitialPointCount-1))

    def evaluate(self, non_deform_coord):
        return pymbolic.evaluate(self.Expression, {"t": non_deform_coord })

    def evaluate_to_vector(self, non_deform_coord):
        result = num.zeros((2,), num.Float)
        result[self.DeformationCoordinate] = pymbolic.evaluate(self.Expression, {"t": non_deform_coord })
        result[1-self.DeformationCoordinate] = non_deform_coord
        return result

    def contains_point(self, point, relative_threshold = 1e-10):
        a,b = self.Interval
        non_dc_value = point[1-self.DeformationCoordinate]
        return min(a,b) <= non_dc_value <= max(a,b) and \
               abs(self.evaluate(non_dc_value) -
                   point[self.DeformationCoordinate]) < relative_threshold * abs(b-a)




class ShapeSection:
    """Describes a closed polygon."""
    def __init__(self, shape_guide_list, tracking_id):
        self.ShapeGuideList = shape_guide_list
        self.TrackingId = tracking_id

    def contains_point(self, point, relative_threshold = 1e-10):
        my_shape_guide_list = self.ShapeGuideList
        # if the last element is just a point, explicitly close the
        # polygon
        if not isinstance(self.ShapeGuideList[-1], ShapeGuide):
            my_shape_guide_list = self.ShapeGuideList[:] + \
                                                        [self.ShapeGuideList[0]]

        last_point = None
        for i in my_shape_guide_list:
            if isinstance(i, ShapeGuide):
                a,b = i.Interval
                start_point = i.evaluate_to_vector(a)
                if last_point is not None:
                    dtl, alpha = toybox.distance_to_line(last_point, start_point-last_point, point)
                    if -relative_threshold <= alpha <= 1+relative_threshold and \
                              dtl <= comp.norm_2(start_point-last_point) * relative_threshold:
                        return True
                if i.contains_point(point, relative_threshold):
                    return True
                last_point = i.end_point()
            else:
                if last_point is not None:
                    dtl, alpha = toybox.distance_to_line(last_point, i-last_point, point)
                    if -relative_threshold <= alpha <= 1+relative_threshold and \
                              dtl <= comp.norm_2(i-last_point) * relative_threshold:
                        return True
                last_point = i
        return False




# Abstract mesh structures ---------------------------------------------------------
class MeshChange:
    def __init__(self, before, after):
        self.MeshBefore = before
        self.MeshAfter = after

    def mesh_before(self):
        return self.MeshBefore

    def mesh_after(self):
        """Returns the changed mesh. Initially, this mesh will not be in
        the generated state.
        """
        return self.MeshAfter

    def change_solution(self, before_solution_vector):
        """Interpolates a solution vector for the "before" mesh into one for
        the "after" mesh, which it returns. The "after" mesh is assumed to be
        in the generated state at the point of the call.
        """
        dofm = self.MeshAfter.dofManager()
        after = num.zeros((self.MeshAfter.dofManager().countDegreesOfFreedom(),), 
            before_solution_vector.typecode())
        for node in self.MeshAfter.nodes():
            before_el = self.MeshBefore.findElement(node.coordinates())
            if before_el is not None:
                after[dofm.getDegreeOfFreedomNumber(node)] = \
                    before_el.getSolutionFunction(before_solution_vector)(node.coordinates())
        return after





class Mesh:
    """Provides a generic base class for any kind of finite element mesh.
    Includes a cache facility for element finders. Once generated,
    meshes are immutable.
    """

    def __init__(self):
        """Sets up parameters and basic data for the mesh and generates it.
        """
        self.DOFManager = element.DOFManager()
        self.Elements = []

    def dimensions(self):
        """Return the number of dimensions of this mesh."""
        raise NotImplementedError

    def dof_manager(self):
        return self.DOFManager

    def elements(self):
        return self.Elements

    def refine(self, element_needs_refining):
        """Returns a MeshChange that represents a refined mesh
        based on the element_needs_refining function passed as
        the argument.

        Returns None if the mesh is not refinable.
        """
        return None

    def find_element(self, point):
        # double-curvy bend: instance-modifying code.
        self.findElement = spatial_btree.build_element_finder(self.Elements)
        return self.findElement(point)




# One-dimensional mesh --------------------------------------------------
class OneDimensionalMesh(Mesh):
    def __init__(self, a, b, n, 
                 element_factory=element.ElementFactory1D(1),
                 left_tracking_id="left", 
                 right_tracking_id="right"):
        Mesh.__init__(self)
        a_vec = num.array([a], num.Float)
        b_vec = num.array([b], num.Float)
        points = comp.linspace(a_vec, b_vec, n)

        self.DOFManager.register_node(0, points[0], left_tracking_id)
        for i, p in list(enumerate(points))[1:-1]:
            self.DOFManager.register_node(i, p)
        self.DOFManager.register_node(-1, points[-1], right_tracking_id)

        assert len(self.DOFManager) == len(points)

        for i in range(len(self.DOFManager) - 1):
            el = element_factory( 
                [self.DOFManager[i], self.DOFManager[i+1]],
                self.DOFManager)
            self.Elements.append(el)

    def dimensions(self):
        return 1




# Pyangle mesh ----------------------------------------------------------
class InverseDeformation:
    def __init__(self, f, f_prime, matinv, nc0):
        self.F = f
        self.FPrime = f_prime
        self.MatrixInverse = matinv
        self.NC0 = nc0

    def __call__(self, y):
        def newton_func(x):
            return self.F(x) - y
        try:
            linear_inverse = self.MatrixInverse*(y - self.NC0)
            return toybox.find_vector_zero_by_newton(newton_func, 
                                                     self.FPrime, linear_inverse)
        except RuntimeError:
            print "WARNING: Transform inversion failed"
            return linear_inverse




class _PyangleMesh(Mesh):
    """This is an internal class.
    Do not use from the outside.
    """
    def __init__(self, generating_parameters,  
                 shape_sections, 
                 element_factory):
        Mesh.__init__(self)
        self.GeneratingParameters = generating_parameters
        self.ElementFactory = element_factory
        self.ShapeSections = shape_sections
        self._postprocess_triangle_output(generating_parameters)

    def find_shape_guide_number(self, number):
        if number == 0:
            return None
        number -= 1
        for sec in self.ShapeSections:
            if number >= len(sec.ShapeGuideList):
                number -= len(sec.ShapeGuideList)
            else:
                return sec.ShapeGuideList[number]

        raise IndexError, "shape guide number out of range"

    def find_shape_section_by_number(self, number):
        if number == 0:
            return None
        number -= 1
        for sec in self.ShapeSections:
            if number >= len(sec.ShapeGuideList):
                number -= len(sec.ShapeGuideList)
            else:
                return sec

        raise IndexError, "shape guide number out of range"

    def dimensions(self):
        return 2

    def _postprocess_triangle_output(self, out_p):
        pts = out_p.Points
        tris = out_p.Triangles

        # read and reposition nodes -----------------------------------------------
        for no in range(pts.size()):
            marker = out_p.PointMarkers.get(no)
            guide = self.find_shape_guide_number(marker)
            section = self.find_shape_section_by_number(marker)

            tracking_id = None
            if section:
                tracking_id = section.TrackingId

            if isinstance(guide, ShapeGuide):
                c = guide.DeformationCoordinate
                pts.set_sub(no, c, guide.evaluate(pts.get_sub(no, 1-c)))

            coordinates = num.array([pts.get_sub(no, 0), pts.get_sub(no, 1)])
            self.DOFManager.register_node(no, coordinates, tracking_id, section)

        # build elements ----------------------------------------------------------
        elements = []
        for tri in range(tris.size()):
            my_nodes = [self.DOFManager.get_node_by_tag(nd) 
                        for nd in [tris.get_sub(tri, 0), 
                                   tris.get_sub(tri, 1),
                                   tris.get_sub(tri, 2)]]

            possible_guides = []

            for node in my_nodes:
                guide = self.find_shape_guide_number(out_p.PointMarkers.get(node.Tag))

                if isinstance(guide, ShapeGuide) and guide not in possible_guides:
                    my_guided_nodes = []
                    for index, no_ in enumerate(my_nodes):
                        if guide.contains_point(no_.Coordinates):
                            my_guided_nodes.append(no_)
                        else:
                            my_index_unguided = index

                    if len(my_guided_nodes) == 3:
                        raise RuntimeError, "can't guide three points with a single guide"
                    elif len(my_guided_nodes) == 2:
                        possible_guides.append(guide)
                        if len(possible_guides) == 1:
                            index_unguided = my_index_unguided

            if len(possible_guides) > 1:
                print "WARNING: found more than one possible guide"

            if len(possible_guides) == 0 or not possible_guides[0].UseExactElements:
                # none of the edges are guided
                el = self.ElementFactory(my_nodes, self.DOFManager)
                elements.append(el)

            elif len(possible_guides) >= 1:
                # one edge is guided

                # rotate list such that unguided node comes first
                my_nodes = [my_nodes[(index_unguided+i)%3] for i in range(3)]

                guide = possible_guides[0]

                # 2
                # |\
                # 0-1
                #
                # 1<->2 is the guided edge
                e_x = pymbolic.subscript(pymbolic.var("x"), 0)
                e_y = pymbolic.subscript(pymbolic.var("x"), 1)
                e_vec = pymbolic.vector.Vector(e_x, e_y)

                # alpha is the blend factor
                # alpha is zero on 2<->0<->1 and one on 1<->2.
                expr_alpha = 4*e_x*e_y/(1-(e_x-e_y)**2)

                # calculate the forward linear transform 
                node_coords = [node.Coordinates for node in my_nodes]
                nc0 = node_coords[0]
                mat = num.array([ n - nc0 for n in node_coords[1:] ]).T
                matinv = 1/mat

                # assemble expression for forward linear tranform
                expr_linear_transform = nc0 + mat * e_vec

                # assemble expression for forward nonlinear transform
                deform_coord = guide.DeformationCoordinate
                non_deform_coord = 1-guide.DeformationCoordinate

                # map (x-y)-interval [-1,1] to [a,b]
                a = node_coords[1][non_deform_coord]
                b = node_coords[2][non_deform_coord]
                expr_guide_argument = (a-b)/2.*(e_x-e_y) + (a+b)/2.
                expr_guide = pymbolic.substitute(
                    guide.Expression, {pymbolic.var("t"): expr_guide_argument})

                expr_transform = list(expr_linear_transform[:])

                expr_reference_x = (1+e_x - e_y)*0.5
                expr_reference_y = (1-(e_x - e_y))*0.5
                expr_linear_reference = nc0[deform_coord] + \
                                        pymbolic.linear_combination(mat[deform_coord], 
                                                                    [expr_reference_x, expr_reference_y])
                expr_transform[deform_coord] = expr_alpha*(expr_guide-expr_linear_reference) \
                                               + expr_linear_transform[deform_coord]

                # compose inverse nonlinear transform
                func_transform = pymbolic.VectorFunction(expr_transform, [pymbolic.var("x")])
                func_transform_jacobian = pymbolic.MatrixFunction(
                    pymbolic.jacobian(expr_transform, e_vec), [pymbolic.var("x")])

                el = self.ElementFactory.make_distorted(
                    my_nodes, self.DOFManager, 
                    func_transform, 
                    func_transform_jacobian,
                    InverseDeformation(func_transform, func_transform_jacobian,
                                       matinv, nc0))
                elements.append(el)

        self.Elements = elements

    def get_refinement(self, element_needs_refining):
        input_p = self.GeneratingParameters.copy()
        
        input_p.TriangleAreas.setup()
        marked_elements = 0
        for i, el in zip(range(len(self.Elements)), self.Elements):
            if element_needs_refining(el):
                input_p.TriangleAreas.set(i, el.area() * 0.7)
                marked_elements += 1
            else:
                input_p.TriangleAreas.set(i, -1)

        if marked_elements == 0:
            raise RuntimeError, "No elements marked for refinement."
        
        return MeshChange(self, _TwoDimensionalRefinedMesh(
            input_p, self.ShapeSections, self.ElementFactory))

class TwoDimensionalMesh(_PyangleMesh):
    def __init__(self, shape_sections, hole_starts=[], 
                 refinement_func=None,
                 element_factory=element.ElementFactory2DTriangle(2)):
        """Specifies the data for an analytically bounded mesh.
        The shape guide list may consist of a list of guides and
        points. 
        """

        points = []
        point_markers = []
        segments = []
        segment_markers = []

        shape_guide_index = 1
        for section in shape_sections:
            section_start_point_index = len(points)
            for shape in section.ShapeGuideList:
                try:
                    a,b = shape.Interval
                    h = float(b-a)/shape.InitialPointCount
                    c = shape.DeformationCoordinate

                    render_point_count = shape.InitialPointCount
                    if shape.RenderFinalPoint:
                        render_point_count += 1
                    for pt_idx in range(render_point_count):
                        pt = num.zeros((2,), num.Float)
                        if pt_idx == shape.InitialPointCount:
                            pt[1-c] = b
                        else:
                            pt[1-c] = a + pt_idx * h
                        pt[c] = pymbolic.evaluate(shape.Expression, {"t": pt[1-c]})
                        segments.append((len(points), len(points) + 1))
                        points.append(pt)
                        if pt_idx == shape.InitialPointCount:
                            point_markers.append(shape_guide_index)
                            segment_markers.append(shape_guide_index)
                        else:
                            point_markers.append(shape_guide_index)
                            segment_markers.append(shape_guide_index)
                except AttributeError:
                    # Oops. It's not a ShapeGuide, must be a point.
                    segments.append((len(points), len(points) + 1))
                    points.append(shape)
                    point_markers.append(shape_guide_index)
                    segment_markers.append(shape_guide_index)
                shape_guide_index += 1

            # bend end of section back to start point
            segments[-1] = (len(points)-1, section_start_point_index)
        
        input_p = pyangle.TriangulationParameters()
        pyangle.set_points(input_p, points, point_markers)
        pyangle.set_segments(input_p, segments, segment_markers)
        pyangle.set_holes(input_p, hole_starts)

        output_p = pyangle.triangulate(input_p, 
                                       refinement_func=refinement_func)
        
        _PyangleMesh.__init__(self, output_p, shape_sections,
                              element_factory)




class _TwoDimensionalRefinedMesh(_PyangleMesh):
    def __init__(self, input_p, shape_sections, element_factory):
        _PyangleMesh.__init__(self, pyangle.refine(input_p), 
                              shape_sections,
                              element_factory)





# Tools -----------------------------------------------------------------
def get_boundary_edges(mesh):
    result = pytools.DictionaryWithDefault(lambda x: [])
    for el in mesh.elements():
        #print el.transformToReal(nd.Coordinates)
        try:
            unit_coords = dict([(nd,el.transform_to_unit(nd.Coordinates)) for nd in el.nodes()])
        except RuntimeError:
            for nd in el.nodes():
                print "ALL", nd.Coordinates
            raise
        center = pytools.average(unit_coords.values())

        def cent_angle(nd):
            my_unit_coords = unit_coords[nd]-center
            return math.atan2(my_unit_coords[1], my_unit_coords[0])
        # sort nodes in counterclockwise order
        sorted_nodes = el.nodes()[:]
        sorted_nodes.sort(lambda n1, n2: cmp(cent_angle(n1), cent_angle(n2)))
        sorted_nodes.append(sorted_nodes[0])

        last_node = sorted_nodes[0]
        for nd in sorted_nodes[1:]:
            ss1 = last_node.ShapeSection
            ss2 = nd.ShapeSection
            if ss1 and \
                      ss1.contains_point(nd.Coordinates) and \
                      ss1.contains_point(last_node.Coordinates):
                result[ss1.TrackingId].append((last_node, nd))
            elif ss2 and \
                      ss2.contains_point(nd.Coordinates) and \
                      ss2.contains_point(last_node.Coordinates):
                result[ss2.TrackingId].append((last_node, nd))
            last_node = nd
    return result


