import pylinear.matrices as num
import element
import matrix_builder
import pyangle
import fempy.spatial_btree as spatial_btree




def _makeDistortedLinearElement(nodes, dof_manager):
  node_coords = [ node.coordinates() for node in nodes ]

  mat = num.array([ n - node_coords[0] for n in node_coords[1:] ])
  matinv = la.inverse(mat)
  nc0 = node_coords[0]

  vars = [ ("variable","0"), ("variable","1") ]
  inv_vars = [ 
    ("-", ("variable","0"), nc0[0]), 
    ("-", ("variable","1"), nc0[1]) ]

  return element.tDistortedTwoDimensionalLinearTriangularFiniteElement(
      nodes, [ 
        ("+", nc0[0], expression.linearCombination(mat[0], vars)), 
        ("+", nc0[1], expression.linearCombination(mat[1], vars)), 
      ], [
        expression.linearCombination(matinv[0], inv_vars), 
        expression.linearCombination(matinv[1], inv_vars), 
      ], dof_manager)




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




class tRectangularMesh(tMesh):
  def __init__(self, nx, ny, x = 1., y = 1., order = 1):
    tMesh.__init__(self)
    self.Steps = [nx, ny]
    self.Extent = [x,y]
    self.Order = 1

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




class tPyangleGeneratedMesh(tMesh):
  def _postprocessTriangleOutput(self, out_p):
    pts = out_p.Points
    tris = out_p.Triangles

    nodes = []
    for node in range(pts.size()):
      nodes.append(element.tNode(num.array([ pts.getSub(node, 0), pts.getSub(node, 1) ])))

    between_nodes = {}

    def between(node1, node2):
      if (node1,node2) in between_nodes:
        return between_nodes[ node1,node2 ]
      else:
        new_node = element.tNode((node1.coordinates() + node2.coordinates()) / 2)
        nodes.append(new_node)
        between_nodes[ node1,node2 ] = \
          between_nodes[ node2,node1 ] = new_node
        return new_node

    elements = []
    for tri in range(tris.size()):
      a = nodes[ tris.getSub(tri, 0) ]
      b = nodes[ tris.getSub(tri, 1) ]
      c = nodes[ tris.getSub(tri, 2) ]
      if self.Order == 2:
        elements.append(element.tTwoDimensionalQuadraticTriangularFiniteElement(
          [ a, b, c, between(a, b), between(b, c), between(c,a) ], 
          self.DOFManager))
      elif self.Order == 1:
        elements.append(element.tTwoDimensionalLinearTriangularFiniteElement([a,b,c], self.DOFManager))
      else:
        raise RuntimeError, "invalid order specified"

    self.Nodes = nodes
    self.Elements = elements

    # find boundary nodes -----------------------------------------------------
    self.BoundaryNodes = []
    for i in range(out_p.PointMarkers.size()):
      if out_p.PointMarkers.get(i):
        self.BoundaryNodes.append(nodes[i])

    self.LastOutputParameters = out_p

  def getRefinement(self, element_needs_refining):
    in_p = self.LastOutputParameters.copy()
    
    for i, el in zip(range(len(self.Elements)), self.Elements):
      if element_needs_refining(el):
        in_p.TriangleAreas.set(i, el.area() * 0.3333)

    print "FIXME: have a look a the area constraint array"
    raw_input()

    return tMeshChange(self, _tTwoDimensionalShapedRefinedMesh(self.Order, in_p))




class tTwoDimensionalShapedMesh(tPyangleGeneratedMesh):
  def __init__(self, shape_points, order = 1, refinement_func = None):
    tMesh.__init__(self)
    self.ShapePoints = shape_points
    self.Order = order
    self.RefinementFunction = refinement_func

  def generate(self):
    if self.Nodes:
      raise RuntimeError, "generate() may not be called twice"

    out_p = pyangle.triangulateArea(self.ShapePoints, 
      refinement_func = self.RefinementFunction)

    # The refinement function may carry expensive references to other
    # meshes. We should not depend on its presence any longer 
    # than necessary.
    del self.RefinementFunction
    
    self._postprocessTriangleOutput(out_p)




class _tTwoDimensionalShapedRefinedMesh(tPyangleGeneratedMesh):
  def __init__(self, order, in_p):
    tMesh.__init__(self)
    self.Order = order
    self.InputParameters = in_p

  def generate(self):
    if self.Nodes:
      raise RuntimeError, "generate() may not be called twice"

    out_p = pyangle.refine(self.InputParameters)
    del self.InputParameters

    self._postprocessTriangleOutput(out_p)





class tAnalyticallyBoundedMesh(tMesh):
  pass
