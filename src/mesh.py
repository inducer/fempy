import pylinear.matrices as num
import pylinear.linear_algebra as la
import element
import expression
import matrix_builder
import pyangle
import tools
import math
import fempy.spatial_btree as spatial_btree




def _makeDistortedLinearElement(nodes, dof_manager):
  node_coords = [ node.coordinates() for node in nodes ]

  mat = num.transpose(num.array([ n - node_coords[0] for n in node_coords[1:] ]))
  matinv = la.inverse(mat)
  nc0 = node_coords[0]

  variables = [ ("variable","0"), ("variable","1") ]

  def inv_func(point):
    return num.matrixmultiply(matinv, point - nc0)
    
  return element.tDistortedTwoDimensionalTriangularFiniteElement(
      nodes, [ 
        ("+", nc0[0], expression.linearCombination(mat[0], variables)), 
        ("+", nc0[1], expression.linearCombination(mat[1], variables)), 
      ], inv_func, dof_manager)




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
    finder = self.MeshBefore.findElement
    dofm = self.MeshAfter.dofManager()
    after = num.zeros((self.MeshAfter.dofManager().countDegreesOfFreedom(),), 
      before_solution_vector.typecode())
    for node in self.MeshAfter.nodes():
      before_el = finder(node.coordinates())
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
    """Sets up parameters and basic data for the mesh, but should *not*
    generate it."""
    self.ElementFinder = None
    self.DOFManager = element.tDOFManager()
    self.Elements = []

  def dimensions(self):
    """Return the number of dimensions of this mesh. Most often, this will be
    either two or three.
    """
    raise RuntimeError, "Not implemented"

  def dofManager(self):
    return self.DOFManager

  def generate(self):
    """Generate the mesh from the parameters given to __init__().
    This method may only be called once on a given instance.
    """
    pass

  def elements(self):
    return self.Elements

  def getRefinement(self, element_needs_refining):
    """Returns a tMeshChange that represents a refined mesh
    based on the element_needs_refining function passed as
    the argument.

    Returns None if the mesh is not refinable.
    """
    return None

  def findElement(self, element):
    # double-curvy bend: instance-modifying code.
    self.findElement = spatial_btree.buildElementFinder(self.Elements)
    return self.findElement(element)




# Pyangle mesh ----------------------------------------------------------
class tShapeGuide:
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

  def containsPoint(self, point, relative_threshold):
    a,b = self.Interval
    return abs(self.evaluate(point[1-self.DeformationCoordinate]) -
               point[self.DeformationCoordinate]) < relative_threshold * abs(b-a)




class tShapeSection:
  """Describes a closed polygon."""
  def __init__(self, shape_guide_list, constraint_id):
    self.ShapeGuideList = shape_guide_list
    self.ConstraintId = constraint_id

  def containsPoint(self, point, relative_threshold = 1e-10):
    my_shape_guide_list = self.ShapeGuideList
    # if the last element is just a point, explicitly close the
    # polygon
    if not isinstance(self.ShapeGuideList[-1], tShapeGuide):
      my_shape_guide_list = self.ShapeGuideList[:] + \
                            [self.ShapeGuideList[0]]

    last_point = my_shape_guide_list[0]
    for i in my_shape_guide_list[1:]:
      if isinstance(i, tShapeGuide):
        if i.containsPoint(point, relative_threshold):
          return True
        a,b = i.Interval
        last_point = i.evaluate(b)
      else:
        dtl, alpha = tools.distanceToLine(last_point, i-last_point, point)
        if -relative_threshold <= alpha <= 1+relative_threshold and \
           dtl <= tools.norm2(i-last_point) * relative_threshold:
          return True
        last_point = i




class _tPyangleMesh(tMesh):
  """This is an internal class.
  Do not use from the outside.
  """
  def __init__(self, input_p, order, shape_sections):
    tMesh.__init__(self)
    self.Order = order
    self.InputParameters = input_p
    self.ShapeSections = shape_sections

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

  def _postprocessTriangleOutput(self, out_p):
    pts = out_p.Points
    tris = out_p.Triangles


    # read and reposition nodes -----------------------------------------------
    for no in range(pts.size()):
      marker = out_p.PointMarkers.get(no)
      guide = self.findShapeGuideNumber(marker)
      section = self.findShapeSectionByNumber(marker)

      constraint_id = None
      if section:
        constraint_id = section.ConstraintId

      if isinstance(guide, tShapeGuide):
        c = guide.DeformationCoordinate
        pts.setSub(no, c, guide.evaluate(pts.getSub(no, 1-c)))

      self.DOFManager.registerNode(no,
                                   num.array([pts.getSub(no, 0), pts.getSub(no, 1)]),
                                   constraint_id, section)

    pyangle.writeGnuplotMesh(",,mesh.data", out_p)

    # build elements ----------------------------------------------------------
    elements = []
    for tri in range(tris.size()):
      my_nodes = [self.DOFManager.getNodeByTag(nd) 
                  for nd in [tris.getSub(tri, 0), 
                             tris.getSub(tri, 1),
                             tris.getSub(tri, 2)]]

      def guideContainsNode(guide, node):
        coords = node.Coordinates
        a,b = guide.Interval
        scale = math.fabs(b-a)
        c = guide.DeformationCoordinate
        if min(a,b) <= coords[1-c] <= max(a,b):
          return math.fabs(guide.evaluate(coords[1-c])-coords[c]) < 1e-10 * scale
        else:
          return False

      possible_guides = []
      guided_nodes = []
      for node in my_nodes:
        guide = self.findShapeGuideNumber(out_p.PointMarkers.get(node.Tag))
        if isinstance(guide, tShapeGuide) and guide not in possible_guides:
          my_guided_nodes = []
          for index, no_ in zip(range(3), my_nodes):
            if guideContainsNode(guide, no_):
              my_guided_nodes.append(no_)
            else:
              my_index_unguided = index

          if len(my_guided_nodes) == 3:
            raise RuntimeError, "can't guide three points with a single guide"
          elif len(my_guided_nodes) == 2:
            possible_guides.append(guide)
            if len(possible_guides) == 1:
              guided_nodes = my_guided_nodes
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
        # t = x+y
        expr_t = ("+", ("variable", "0"), ("variable", "1"))
        # s = x-y
        expr_s = ("-", ("variable", "0"), ("variable", "1"))

        # alpha is the blend factor
        # alpha is zero on 2<->0<->1 and one on 1<->2.
        # alpha = (s^2-t^2)/(s^2-1)
        expr_alpha = ("/", ("-", ("**", expr_s, 2), ("**", expr_t, 2)),
                     ("+", 1e-50, ("-", ("**", expr_s, 2), 1)))

        # calculate the forward linear transform 
        node_coords = [node.Coordinates for node in my_nodes]
        nc0 = node_coords[0]
        mat = num.transpose(num.array([ n - nc0 for n in node_coords[1:] ]))
        matinv = la.inverse(mat)

        variables = [ ("variable","0"), ("variable","1") ]

        # assemble expression for forward linear tranform
        expr_linear_transform = [
          ("+", nc0[0], expression.linearCombination(mat[0], variables)),
          ("+", nc0[1], expression.linearCombination(mat[1], variables))]

        # assemble expression for forward nonlinear transform
        deform_coord = guide.DeformationCoordinate
        non_deform_coord = 1-guide.DeformationCoordinate

        # map s-space [-1,1] to [a,b]
        a = node_coords[1][non_deform_coord]
        b = node_coords[2][non_deform_coord]
        expr_guide_argument = ("+",("*",(a-b)/2., expr_s), (b+a)/2.)
        expr_guide = expression.substitute(guide.Expression, 
                                            {"t": expr_guide_argument})
        expr_transform = expr_linear_transform[:]

        expr_x_minus_y = ("-",("variable","0"),("variable","1"))
        expr_reference_x = ("*", ("+", expr_x_minus_y, 1), 0.5)
        expr_reference_y = ("*", ("-", 1, expr_x_minus_y), 0.5)
        expr_linear_reference = ("+", 
                                 nc0[deform_coord], 
                                 expression.linearCombination(mat[deform_coord], 
                                                              [expr_reference_x, expr_reference_y]))
        expr_transform[deform_coord] = ("+",
                                        ("*", expr_alpha, 
                                         ("-", expr_guide, expr_linear_reference)),
                                        expr_linear_transform[deform_coord])

        if False:
          txtoreal = expression.compileVectorField(expr_transform)
          alpha_c = expression.compileScalarField(expr_alpha)
          offset_c = expression.compileScalarField(("-", expr_guide, expr_linear_transform[deform_coord]))

          deform_file = file(",,deform.data", "w")
          segments = 20
          h = 1./segments
          for y_i in range(segments):
            y = float(y_i)/segments
            for x_i in range(segments-y_i):
              x = float(x_i)/segments

              def write(x,y):
                alpha = alpha_c([x,y])
                offset = offset_c([x,y])
                deform = txtoreal([x,y])
                deform_file.write("%f\t%f\t%f\n" % (deform[0], deform[1], alpha))

              # 0
              # 12
              write(x,y+h)
              write(x,y)
              write(x+h,y)
              write(x,y+h)
              deform_file.write("\n\n")

          nodes_before = [[0,0], [1,0], [0,1]]
          for i in range(3):
            assert tools.norm2(node_coords[i] - txtoreal(nodes_before[i])) < 1e-10

        # compose inverse nonlinear transform
        func_transform = expression.compileVectorField(expr_transform)
        func_transform_jacobian = expression.compileMatrixFunction(
          expression.jacobian(expr_transform, ["0", "1"]))
        
        def inv_func(point):
          def newton_func(x):
            return func_transform(x)-point
          return tools.findVectorZeroByNewton(newton_func, 
                                              func_transform_jacobian, 
                                              num.matrixmultiply(matinv, point - nc0))

        elements.append(
          element.tDistortedTwoDimensionalTriangularFiniteElement(
          my_nodes, self.DOFManager, element.QuadraticFormFunctionKit,
          func_transform, 
          func_transform_jacobian,
          inv_func))

    self.Elements = elements
    self.LastOutputParameters = out_p

  def getRefinement(self, element_needs_refining):
    input_p = self.LastOutputParameters.copy()
    
    input_p.TriangleAreas.setup()
    marked_elements = 0
    for i, el in zip(range(len(self.Elements)), self.Elements):
      if element_needs_refining(el):
        input_p.TriangleAreas.set(i, el.area() * 0.3)
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

    points_file = file(",,points.data", "w")

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
            points_file.write("%f\t%f\n" % (pt[0], pt[1]))
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
    _tPyangleMesh.__init__(self, input_p, order, shape_sections)

    self.RefinementFunction = refinement_func

  def generate(self):
    if self.Elements:
      raise RuntimeError, "generate() may not be called twice"

    self._postprocessTriangleOutput(
      pyangle.triangulate(self.InputParameters, refinement_func = self.RefinementFunction))

    # The refinement function may carry expensive references to other
    # meshes. We should not depend on its presence any longer 
    # than necessary.
    del self.RefinementFunction




class _tTwoDimensionalRefinedMesh(_tPyangleMesh):
  def __init__(self, input_p, order, shape_sections):
    _tPyangleMesh.__init__(self, input_p, order, shape_sections)

  def generate(self):
    if self.Elements:
      raise RuntimeError, "generate() may not be called twice"

    self._postprocessTriangleOutput(pyangle.refine(self.InputParameters))
    del self.InputParameters





