import pylinear.matrices as num
import element
import matrix_builder




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
    # FIXME: implement me using the element finders
    raise RuntimeError, "not yet implemented"




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

  def dofManager(self):
    return self.DOFManager

  def generate(self):
    """Generate the mesh from the parameters given to __init__().
    This method may only be called once on a given instance.
    """
    pass

  def getData(self):
    """Returns a tuple (nodes, elements, hanging_nodes) where

    nodes is a list of tNodes

    elements is a list of tElements

    hanging_nodes is a list of three-tuples that each look like
    (node1, hanging_node, node2), where hanging_node's value is fixed
    halfway between node1 and node2.
    """
    return self.Nodes, self.Elements, self.HangingNodes

  def getRefinement(self, mark_for_refinement):
    """Returns a tMeshChange that takes into account the 
    knowledge contained in the error_estimator and this mesh. 
    The error_estimator is a function that receives
    an element as the input and returns some real-valued positive measure 
    of the "error" that is to be eliminated by refinement.

    Returns None if the mesh is not refinable.
    """
    return None

  def getElementFinder(self):
    if self.ElementFinder is None:
      self.ElementFinder = spatial_btree.buildElementFinder(self.Elements)
    return self.ElementFinder




class tRectangularMesh(tMesh):
  def __init__(self, nx, ny, x = 1., y = 1., order = 1):
    tMesh.__init__(self)
    self.Steps = [nx, ny]
    self.Extent = [x,y]
    self.Order = 1
    self.generate()

  def generate(self):
    nx = self.Steps[0]
    ny = self.Steps[1]
    x = self.Extent[0]
    y = self.Extent[1]
    dx = x / nx
    dy = y / ny

    # build nodes
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

    # build elements, pay attention to mathematically positive orientation
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




class tShapedMesh(tMesh):
  def __init__(self, shape_points, order = 1, refinement_func = None):
    tMesh.__init__(self)
    self.ShapePoints = shape_points
    self.Order = order

    def default_refinement_func(v1, v2, v3, area):
      return False

    self.RefinementFunction = refinement_func

  def generate(self):
    if self.Nodes:
      return

    self.DOFManager = matrix_builder.tDOFManager()
    out_p = pyangle.triangulateArea(shape_points, self.RefinementFunction)

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
      if second_order:
        elements.append(element.tTwoDimensionalQuadraticTriangularFiniteElement(
          [ a, b, c, between(a, b), between(b, c), between(c,a) ], 
          self.DOFManager))
      else:
        elements.append(element.tTwoDimensionalLinearTriangularFiniteElement([a,b,c], self.DOFManager))

    self.Nodes = nodes
    self.Elements = elements

    # The refinement function may carry expensive references to other
    # meshes. We should not depend on its presence any longer 
    # than necessary.
    del self.RefinementFunction




class tAnalyticallyBoundedMesh(tMesh):
  pass
