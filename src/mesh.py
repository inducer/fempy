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
    
  return element.tDistortedTwoDimensionalLinearTriangularFiniteElement(
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
    self.DOFManager = matrix_builder.tDOFManager()

    self.Nodes = []
    self.Elements = []
    self.HangingNodes = []
    self.BoundaryNodes = []

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

  def nodes(self):
    return self.Nodes

  def elements(self):
    return self.Elements

  def hangingNodes(self):
    """Returns a list of three-tuples that each look like
    (node1, hanging_node, node2), where hanging_node's value should 
    be fixed halfway between node1 and node2.
    """
    return self.HangingNodes

  def boundaryNodes(self):
    """Returns a list of nodes that form the outer boundary of the given
    geometry.
    """
    return self.BoundaryNodes

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




# rectangular mesh ------------------------------------------------------------
class tRectangularMesh(tMesh):
  def __init__(self, nx, ny, x = 1., y = 1., order = 1):
    tMesh.__init__(self)
    self.Steps = [nx, ny]
    self.Extent = [x,y]
    self.Order = 1

  def dimensions(self):
    return 2

  def generate(self):
    nx = self.Steps[0]
    ny = self.Steps[1]
    x = self.Extent[0]
    y = self.Extent[1]
    dx = x / nx
    dy = y / ny

    # build nodes -------------------------------------------------------------
    all_nodes = []
    nodes = []
    for node_y in range(0, ny + 1):
      line_nodes = []
      for node_x in range(0, nx + 1):
        line_nodes.append(element.tNode(num.array([ node_x * dx, node_y * dy ])))
      nodes.append(line_nodes)
      all_nodes += line_nodes

    between_nodes = {}

    def between(node1, node2):
      if (node1,node2) in between_nodes:
        return between_nodes[ node1,node2 ]
      else:
        new_node = element.tNode((node1.coordinates() + node2.coordinates()) / 2)
        all_nodes.append(new_node)
        between_nodes[ node1,node2 ] = \
          between_nodes[ node2,node1 ] = new_node
        return new_node

    # build elements, pay attention to mathematically positive orientation ----
    elements = []
    for el_y in range(0, ny):
      for el_x in range(0, nx):
        # d c
        # a b
        a = nodes[el_y][el_x]
        b = nodes[el_y][el_x + 1]
        c = nodes[el_y + 1][el_x + 1]
        d = nodes[el_y + 1][el_x]

        if self.Order == 2:
          lower_el = element.tTwoDimensionalQuadraticTriangularFiniteElement(
            [ a,b,d, between(a, b), between(b, d), between(d, a) ], 
            self.DOFManager)
          upper_el = element.tTwoDimensionalQuadraticTriangularFiniteElement(
            [ c,d,b, between(c, d), between(d, b), between(b, c) ], 
            self.DOFManager)
        else:
          lower_el = element.tTwoDimensionalLinearTriangularFiniteElement([ a,b,d ], self.DOFManager)
          upper_el = element.tTwoDimensionalLinearTriangularFiniteElement([ c,d,b ], self.DOFManager)

        elements.append(lower_el)
        elements.append(upper_el)
        
    self.Nodes = all_nodes
    self.Elements = elements

    # find boundary nodes -----------------------------------------------------
    # FIXME: Hello, inefficiency!
    def isEdgeNode(node):
      x = node.coordinates()
      for node_c in node.coordinates():
        if node_c == 0:
          return True
      for node_c, bound_c in node.coordinates(), self.Extent:
        if node_c == bound_c:
          return True
      
    self.BoundaryNodes = filter(isEdgeNode, new_mesh.nodes())




# Pyangle mesh ----------------------------------------------------------
class tShapeGuide:
  def __init__(self, deformation_coordinate, interval, expression, 
      initial_point_count = 3, render_final_point = False):
    self.DeformationCoordinate = deformation_coordinate
    self.Interval = interval
    self.Expression = expression
    self.InitialPointCount = initial_point_count
    self.RenderFinalPoint = render_final_point

class tShapeSection:
  def __init__(self, shape_guide_list, is_boundary):
    self.ShapeGuideList = shape_guide_list
    self.IsBoundary = is_boundary




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
    nodes = []
    boundary_nodes = []

    boundary_file = file(",,boundary.data", "w")
    
    for no in range(pts.size()):
      marker = out_p.PointMarkers.get(no)
      guide = self.findShapeGuideNumber(marker)
      section = self.findShapeSectionByNumber(marker)

      if isinstance(guide, tShapeGuide):
        c = guide.DeformationCoordinate
        pts.setSub(no, c, 
                   expression.evaluate(guide.Expression, {"t": pts.getSub(no, 1-c) }))

      node = element.tNode(num.array([pts.getSub(no, 0), pts.getSub(no, 1)]))

      nodes.append(node)

      if section and section.IsBoundary:
        boundary_nodes.append(node)
        boundary_file.write("%f\t%f\n" % (node.coordinates()[0], node.coordinates()[1]))

    pyangle.writeGnuplotMesh(",,mesh.data", out_p)

    # build elements ----------------------------------------------------------
    elements = []
    for tri in range(tris.size()):
      node_numbers = [tris.getSub(tri, 0), tris.getSub(tri, 1),
        tris.getSub(tri, 2)]
      my_nodes = [nodes[no] for no in node_numbers]

      shape_guides = []
      actual_shape_guides = []
      guided_nodes = []

      for index,no in zip(range(3),node_numbers):
        guide = self.findShapeGuideNumber(out_p.PointMarkers.get(no))
        shape_guides.append(guide)
        if isinstance(guide, tShapeGuide):
          actual_shape_guides.append(guide)
          guided_nodes.append(no)
        else:
          index_unguided = index

      if len(guided_nodes) <= 1 or not tools.allEqual(actual_shape_guides):
        # none of the edges are guided
        elements.append(
            element.tTwoDimensionalLinearTriangularFiniteElement(
              my_nodes, self.DOFManager))
      elif len(guided_nodes) == 2:
        # one edge is guided

        # rotate list such that unguided node comes first
        node_numbers = [node_numbers[(index_unguided+i)%3] for i in range(3)]

        guide = actual_shape_guides[0]

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

        # alpha is zero on 2<->0<->1 and one on 1<->2.
        # alpha = (s^2-t^2)/(s^2-1)
        expr_alpha = ("/", ("-", ("**", expr_s, 2), ("**", expr_t, 2)),
                     ("-", ("**", expr_s, 2), 1))

        # calculate the forward linear transform 
        node_coords = [nodes[no].coordinates() for no in node_numbers]
        mat = num.transpose(num.array([ n - node_coords[0] for n in node_coords[1:] ]))
        matinv = la.inverse(mat)
        nc0 = node_coords[0]

        variables = [ ("variable","0"), ("variable","1") ]

        # assemble expression for forward linear tranform
        expr_linear_transform = [
          ("+", nc0[0], expression.linearCombination(mat[0], variables)),
          ("+", nc0[1], expression.linearCombination(mat[1], variables))]

        # assemble expression for forward nonlinear transform
        deform_coord = guide.DeformationCoordinate
        non_deform_coord = guide.DeformationCoordinate

        expr_guide = expression.substitute(guide.Expression, 
                                            {"t": expr_linear_transform[non_deform_coord]})
        expr_transform = expr_linear_transform[:]
        expr_transform[deform_coord] = ("+",("*",("-", 1, expr_alpha), 
                                             expr_linear_transform[deform_coord]),
                                        ("*", expr_alpha, expr_guide))

        # compose inverse nonlinear transform
        deformed_coordinate_function = expression.compileScalarField(expr_transform[deform_coord])
        deformed_coordinate_function_prime = expression.compileScalarField(
          expression.simplify(
          expression.differentiate(expr_transform[deform_coord], "%d" % deform_coord)))
        non_deform_matrix_inverse = matinv[non_deform_coord]
        
        def inv_func(point):
          non_deform_value = num.innerproduct(non_deform_matrix_inverse, point)
          target_value = point[deform_coord]
          def newton_func(deform_value):
            my_point = num.zeros((2,), num.Float)
            my_point[deform_coord] = deform_value
            my_point[non_deform_coord] = non_deform_value
            return deformed_coordinate_function(my_point) - target_value
          def newton_func_prime(deform_value):
            my_point = num.zeros((2,), num.Float)
            my_point[deform_coord] = deform_value
            my_point[non_deform_coord] = non_deform_value
            return deformed_coordinate_function_prime(my_point)
          return tools.findZeroByNewton(newton_func, newton_func_prime, point[deform_coord])

        elements.append(
          element.tDistortedTwoDimensionalLinearTriangularFiniteElement(
          my_nodes, expr_transform, inv_func, self.DOFManager))
      elif len(guided_nodes) == 3:
        # more than one edge is guided
        raise RuntimeError, "More than one edge is guided. This is currently unsupported."

    self.Nodes = nodes
    self.Elements = elements
    self.BoundaryNodes = boundary_nodes
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
    if self.Nodes:
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
    if self.Nodes:
      raise RuntimeError, "generate() may not be called twice"

    self._postprocessTriangleOutput(pyangle.refine(self.InputParameters))
    del self.InputParameters





