import pylinear.matrices as num
import pylinear.linear_algebra as la
import pylinear.matrix_tools as mtools
import element
import expression
import pyangle
import tools
import math
import fempy.spatial_btree as spatial_btree
import expression_operators as eo




# Mesh description data structures ------------------------------------------------
class tShapeGuide:
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

  def evaluate(self, non_deform_coord):
    return expression.evaluate(self.Expression, {"t": non_deform_coord })

  def evaluateToVector(self, non_deform_coord):
    result = num.zeros((2,), num.Float)
    result[self.DeformationCoordinate] = expression.evaluate(self.Expression, {"t": non_deform_coord })
    result[1-self.DeformationCoordinate] = non_deform_coord
    return result

  def containsPoint(self, point, relative_threshold = 1e-10):
    a,b = self.Interval
    non_dc_value = point[1-self.DeformationCoordinate]
    return min(a,b) <= non_dc_value <= max(a,b) and \
           abs(self.evaluate(non_dc_value) -
               point[self.DeformationCoordinate]) < relative_threshold * abs(b-a)




class tShapeSection:
  """Describes a closed polygon."""
  def __init__(self, shape_guide_list, tracking_id):
    self.ShapeGuideList = shape_guide_list
    self.TrackingId = tracking_id

  def containsPoint(self, point, relative_threshold = 1e-10):
    my_shape_guide_list = self.ShapeGuideList
    # if the last element is just a point, explicitly close the
    # polygon
    if not isinstance(self.ShapeGuideList[-1], tShapeGuide):
      my_shape_guide_list = self.ShapeGuideList[:] + \
                            [self.ShapeGuideList[0]]

    last_point = None
    for i in my_shape_guide_list:
      if isinstance(i, tShapeGuide):
        a,b = i.Interval
        start_point = i.evaluateToVector(a)
        if last_point is not None:
          dtl, alpha = tools.distanceToLine(last_point, start_point-last_point, point)
          if -relative_threshold <= alpha <= 1+relative_threshold and \
               dtl <= tools.norm2(start_point-last_point) * relative_threshold:
            return True
        if i.containsPoint(point, relative_threshold):
          return True
        a,b = i.Interval
        last_point = i.evaluateToVector(b)
      else:
        if last_point is not None:
          dtl, alpha = tools.distanceToLine(last_point, i-last_point, point)
          if -relative_threshold <= alpha <= 1+relative_threshold and \
               dtl <= tools.norm2(i-last_point) * relative_threshold:
            return True
        last_point = i
    return False




# Abstract mesh structures ---------------------------------------------------------
class tMeshChange:
  def __init__(self, before, after):
    self.MeshBefore = before
    self.MeshAfter = after

  def meshBefore(self):
    return self.MeshBefore

  def meshAfter(self):
    """Returns the changed mesh. Initially, this mesh will not be in
    the generated state.
    """
    return self.MeshAfter

  def changeSolution(self, before_solution_vector):
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





class tMesh:
  """Provides a generic base class for any kind of finite element mesh.
  Includes a cache facility for element finders. Once generated,
  meshes are immutable.
  """

  def __init__(self):
    """Sets up parameters and basic data for the mesh and generates it.
    """
    self.ElementFinder = None
    self.DOFManager = element.tDOFManager()
    self.Elements = []

  def dimensions(self):
    """Return the number of dimensions of this mesh."""
    raise NotImplementedError

  def dofManager(self):
    return self.DOFManager

  def elements(self):
    return self.Elements

  def getRefinement(self, element_needs_refining):
    """Returns a tMeshChange that represents a refined mesh
    based on the element_needs_refining function passed as
    the argument.

    Returns None if the mesh is not refinable.
    """
    return None

  def findElement(self, point):
    # double-curvy bend: instance-modifying code.
    self.findElement = spatial_btree.buildElementFinder(self.Elements)
    return self.findElement(point)




# One-dimensional mesh --------------------------------------------------
class tOneDimensionalMesh(tMesh):
  def __init__(self, a, b, n, left_tracking_id, right_tracking_id):
    tMesh.__init__(self)
    a_vec = num.array([a], num.Float)
    b_vec = num.array([b], num.Float)
    points = mtools.linspace(a_vec, b_vec, n)

    self.DOFManager.registerNode(0, points[0], left_tracking_id)
    for i, p in list(enumerate(points))[1:-1]:
      self.DOFManager.registerNode(i, p)
    self.DOFManager.registerNode(-1, points[-1], right_tracking_id)

    assert len(self.DOFManager) == len(points)

    for i in range(len(self.DOFManager) - 1):
      self.Elements.append(element.tOneDimensionalFiniteElement(
        [self.DOFManager[i], self.DOFManager[i+1]],
        self.DOFManager,
        element.LinearFormFunctionKit1D))

  def dimensions(self):
    return 1




# Pyangle mesh ----------------------------------------------------------
class tInverseDeformation:
  def __init__(self, f, f_prime, matinv, nc0):
    self.F = f
    self.FPrime = f_prime
    self.MatrixInverse = matinv
    self.NC0 = nc0

  def __call__(self, y):
    def newton_func(x):
      return self.F(x) - y
    try:
      linear_inverse = num.matrixmultiply(self.MatrixInverse, y - self.NC0)
      return tools.findVectorZeroByNewton(newton_func, 
                                          self.FPrime, 
                                          linear_inverse)
    except RuntimeError:
      print "WARNING: Transform inversion failed"
      return linear_inverse




class _tPyangleMesh(tMesh):
  """This is an internal class.
  Do not use from the outside.
  """
  def __init__(self, generating_parameters, order, shape_sections):
    tMesh.__init__(self)
    self.GeneratingParameters = generating_parameters
    self.Order = order
    self.ShapeSections = shape_sections
    self.__postprocessTriangleOutput__(generating_parameters)

  def __getstate__(self):
    my_dict = self.__dict__.copy()
    del my_dict["ElementFinder"]
    return my_dict

  def __setstate__(self, my_dict):
    tMesh.__init__(self)
    self.__dict__.update(my_dict)

  def findShapeGuideNumber(self, number):
    if number == 0:
      return None
    number -= 1
    for sec in self.ShapeSections:
      if number >= len(sec.ShapeGuideList):
        number -= len(sec.ShapeGuideList)
      else:
        return sec.ShapeGuideList[number]

    raise IndexError, "shape guide number out of range"

  def findShapeSectionByNumber(self, number):
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

  def __postprocessTriangleOutput__(self, out_p):
    pts = out_p.Points
    tris = out_p.Triangles

    # read and reposition nodes -----------------------------------------------
    for no in range(pts.size()):
      marker = out_p.PointMarkers.get(no)
      guide = self.findShapeGuideNumber(marker)
      section = self.findShapeSectionByNumber(marker)

      tracking_id = None
      if section:
        tracking_id = section.TrackingId

      if isinstance(guide, tShapeGuide):
        c = guide.DeformationCoordinate
        pts.setSub(no, c, guide.evaluate(pts.getSub(no, 1-c)))

      coordinates = num.array([pts.getSub(no, 0), pts.getSub(no, 1)])
      self.DOFManager.registerNode(no, coordinates, tracking_id, section)

    # build elements ----------------------------------------------------------
    elements = []
    for tri in range(tris.size()):
      my_nodes = [self.DOFManager.getNodeByTag(nd) 
                  for nd in [tris.getSub(tri, 0), 
                             tris.getSub(tri, 1),
                             tris.getSub(tri, 2)]]

      possible_guides = []

      for node in my_nodes:
        guide = self.findShapeGuideNumber(out_p.PointMarkers.get(node.Tag))

        if isinstance(guide, tShapeGuide) and guide not in possible_guides:
          my_guided_nodes = []
          for index, no_ in enumerate(my_nodes):
            if guide.containsPoint(no_.Coordinates):
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
        elements.append(
          element.tTwoDimensionalTriangularFiniteElement(
          my_nodes, self.DOFManager, element.QuadraticFormFunctionKit))
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

        # t and s are substitution variables used within alpha
        e_x = (eo.VARIABLE, "0")
        e_y = (eo.VARIABLE, "1")
        expr_x_minus_y = (eo.MINUS,e_x,e_y)

        # alpha is the blend factor
        # alpha is zero on 2<->0<->1 and one on 1<->2.
        # alpha = -4xy/((x-y)^2-1)
        expr_alpha = (eo.DIVIDE, (eo.NEG, (eo.TIMES, (eo.TIMES, 4, e_x), e_y)), 
                      (eo.PLUS, 1e-50, 
                       (eo.MINUS, (eo.POWER, expr_x_minus_y, 2), 
                        1)))

        # calculate the forward linear transform 
        node_coords = [node.Coordinates for node in my_nodes]
        nc0 = node_coords[0]
        mat = num.transpose(num.array([ n - nc0 for n in node_coords[1:] ]))
        matinv = la.inverse(mat)

        variables = [e_x, e_y]

        # assemble expression for forward linear tranform
        expr_linear_transform = [
          (eo.PLUS, nc0[0], expression.linearCombination(mat[0], variables)),
          (eo.PLUS, nc0[1], expression.linearCombination(mat[1], variables))]

        # assemble expression for forward nonlinear transform
        deform_coord = guide.DeformationCoordinate
        non_deform_coord = 1-guide.DeformationCoordinate

        # map (x-y)-interval [-1,1] to [a,b]
        a = node_coords[1][non_deform_coord]
        b = node_coords[2][non_deform_coord]
        expr_guide_argument = (eo.PLUS,(eo.TIMES,(a-b)/2., expr_x_minus_y), 
                               (a+b)/2.)
        expr_guide = expression.substitute(guide.Expression, 
                                            {"t": expr_guide_argument})
        expr_transform = expr_linear_transform[:]

        expr_reference_x = (eo.TIMES, (eo.PLUS, expr_x_minus_y, 1), 0.5)
        expr_reference_y = (eo.TIMES, (eo.MINUS, 1, expr_x_minus_y), 0.5)
        expr_linear_reference = (eo.PLUS, nc0[deform_coord], 
                                 expression.linearCombination(mat[deform_coord], 
                                                              [expr_reference_x, expr_reference_y]))
        expr_transform[deform_coord] = (eo.PLUS,
                                        (eo.TIMES, expr_alpha, 
                                         (eo.MINUS, expr_guide, 
                                          expr_linear_reference)),
                                        expr_linear_transform[deform_coord])

        # compose inverse nonlinear transform
        func_transform = expression.compileVectorField(expr_transform)
        func_transform_jacobian = expression.compileMatrixFunction(
          expression.jacobian(expr_transform, ["0", "1"]))
        
        elements.append(
          element.tDistortedTwoDimensionalTriangularFiniteElement(
          my_nodes, self.DOFManager, element.QuadraticFormFunctionKit,
          func_transform, 
          func_transform_jacobian,
          tInverseDeformation(func_transform, func_transform_jacobian,
                              matinv, nc0)))

    self.Elements = elements

  def getRefinement(self, element_needs_refining):
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
    
    return tMeshChange(self, _tTwoDimensionalRefinedMesh(
      input_p, self.Order, self.ShapeSections))

class tTwoDimensionalMesh(_tPyangleMesh):
  def __init__(self, shape_sections, hole_starts = [], order = 1, refinement_func = None):
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
            pt[c] = expression.evaluate(shape.Expression, {"t": pt[1-c]})
            segments.append((len(points), len(points) + 1))
            points.append(pt)
            if pt_idx == shape.InitialPointCount:
              point_markers.append(shape_guide_index)
              segment_markers.append(shape_guide_index)
            else:
              point_markers.append(shape_guide_index)
              segment_markers.append(shape_guide_index)
        except AttributeError:
          # Oops. It's not a tShapeGuide, must be a point.
          segments.append((len(points), len(points) + 1))
          points.append(shape)
          point_markers.append(shape_guide_index)
          segment_markers.append(shape_guide_index)
        shape_guide_index += 1

      # bend end of section back to start point
      segments[-1] = (len(points)-1, section_start_point_index)
    
    input_p = pyangle.tTriangulationParameters()
    pyangle.setPoints(input_p, points, point_markers)
    pyangle.setSegments(input_p, segments, segment_markers)
    pyangle.setHoles(input_p, hole_starts)
    _tPyangleMesh.__init__(self, pyangle.triangulate(input_p, refinement_func = refinement_func), order, shape_sections)




class _tTwoDimensionalRefinedMesh(_tPyangleMesh):
  def __init__(self, input_p, order, shape_sections):
    _tPyangleMesh.__init__(self, pyangle.refine(input_p), order, shape_sections)





# Tools -----------------------------------------------------------------
def getBoundaryEdges(mesh):
  result = tools.tDictionaryWithDefault(lambda x: [])
  for el in mesh.elements():
    #print el.transformToReal(nd.Coordinates)
    try:
      unit_coords = dict([(nd,el.transformToUnit(nd.Coordinates)) for nd in el.nodes()])
    except RuntimeError:
      for nd in el.nodes():
        print "ALL", nd.Coordinates
      raise
    center = tools.average(unit_coords.values())

    def cent_angle(nd):
      my_unit_coords = unit_coords[nd]-center
      return math.atan2(my_unit_coords[1], my_unit_coords[0])
    # sort nodes in counterclockwise order
    sorted_nodes = el.nodes()[:]
    sorted_nodes.sort(lambda n1, n2: cmp(cent_angle(n1), cent_angle(n2)))
    sorted_nodes.append(sorted_nodes[0])

    last_node = sorted_nodes[0]
    last_shape_section = sorted_nodes[0].ShapeSection
    for nd in sorted_nodes[1:]:
      ss1 = last_node.ShapeSection
      ss2 = nd.ShapeSection
      if ss1 and \
           ss1.containsPoint(nd.Coordinates) and \
           ss1.containsPoint(last_node.Coordinates):
        result[ss1.TrackingId].append((last_node, nd))
      elif ss2 and \
           ss2.containsPoint(nd.Coordinates) and \
           ss2.containsPoint(last_node.Coordinates):
        result[ss2.TrackingId].append((last_node, nd))
      last_node = nd
      last_shape_section = nd.ShapeSection
  return result


