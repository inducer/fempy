import pylinear.matrices as num




class tVisualizationData:
  """This object represents the visualization for one finite element.
  by specifying a number of triangles and the values of the solution
  function at their nodes. This is what the data mean:

  Nodes is a list of tNodes.

  NodeValues specifies the values of the solution function at the
  nodes given by NodeNumbers, in the same order.

  Triangles is a list of lists of node numbers. To refer to the
  above NodeNumbers array, use regular-zero based indexing into that
  array. But an element may introduce extra nodes into its 
  visualization that are not necessarily mentioned in the DOF manager.
  These are given in the Extra* arrays and referred to by negative
  indices in Triangles. (where -1 -> 0, -2 -> 1 and so on)

  ExtraNodeCoordinates is an array of Numeric arrays of (n,) shape
  (i.e. vectors), giving the coordinates of the extra nodes.

  ExtraNodeValues gives the values of the solution function at the
  extra nodes.

  ExtraPolygons specifies a list of polygons in the same format as
  `Triangles' that allows to specify a 'coarser view' of the given mesh.
  This may be used to show element boundaries. Not all backends show
  this extra information.
  """

  def __init__(self, nodes, node_values, triangles,
      extra_node_coordinates = [], extra_node_values = [], extra_polygons = []):
    self.Nodes = nodes
    self.NodeValues = node_values
    self.Triangles = triangles
    self.ExtraNodeCoordinates = extra_node_coordinates
    self.ExtraNodeValues = extra_node_values
    self.ExtraPolygons = extra_polygons




def _compileInfo(elements, node_data):
  node_number_map = {}
  nodes = []
  values = []
  tris = []
  extra_polys = []

  for el in elements:
    data = el.visualizationData(node_data)
    local_nodes = {}

    index = 0
    for node in data.Nodes:
      node_number = node.Number
      if node_number in node_number_map:
	local_nodes[ index ] = node_number_map[ node_number ]
      else:
	local_nodes[ index ] = len(nodes)
	node_number_map[ node ] = len(nodes)
	nodes.append(node.Coordinates)
	values.append(data.NodeValues[ index ])
      index += 1

    index = -1
    for coords, value in zip(data.ExtraNodeCoordinates, data.ExtraNodeValues):
      local_nodes[ index ] = len(nodes)
      nodes.append(coords)
      values.append(value)
      index -= 1

    for na,nb,nc in data.Triangles:
      tris.append((local_nodes[na], local_nodes[nb], local_nodes[nc]))

    for polygon in data.ExtraPolygons:
      extra_polys.append([local_nodes[node] for node in polygon])

  return nodes, values, tris, extra_polys



def _writeGnuplotFile(name, info):
  gnuplot_file = file(name, "w")

  def writeNode(node):
    gnuplot_file.write("%f %f %f\n" % (
	nodes[node][0],
	nodes[node][1],
	values[node]))

  nodes,values,triangles,extra_polys = info
  for tri_nodes in triangles:
    for node in tri_nodes:
      writeNode(node)
    writeNode(tri_nodes[0])
    gnuplot_file.write("\n\n")




def _writeVtkFile(name, name_grid, info):
  import pyvtk

  nodes,values,triangles,extra_polys = info

  my_nodes = []
  for node in nodes:
    my_nodes.append([ node[0], node[1], 0 ])

  structure = pyvtk.PolyData(points = my_nodes, polygons = triangles)
  structure_grid = pyvtk.PolyData(points = my_nodes, polygons = extra_polys)

  pointdata = pyvtk. PointData(
      pyvtk. Scalars(values, name="node_data", lookup_table = "default"))

  vtk = pyvtk.VtkData(structure, "FEM result triangles", pointdata)
  vtk.tofile(name, "ascii")

  vtk = pyvtk.VtkData(structure_grid, "FEM result grid", pointdata)
  vtk.tofile(name_grid, "ascii")



def _writeMatlabFile(name, info):
  m_file = file(name, "w")

  def writeMatlabVector(name, data):
    (h,) = data.shape
    m_file.write("%s = [\n" % name)
    for i in range(0,h):
      if i != h-1:
	m_file.write("%f;\n" % data[ i ])
      else:
	m_file.write("%f" % data[ i ])
    m_file.write("];\n")

  def writeMatlabMatrix(name, data):
    (h,w) = data.shape
    m_file.write("%s = [\n" % name)
    for i in range(0,h):
      for j in range(0,w):
	if j != w-1:
	  m_file.write("%f," % data[ i,j ])
	else:
	  if i != h-1:
	    m_file.write("%f;\n" % data[ i,j ])
	  else:
	    m_file.write("%f" % data[ i,j ])
    m_file.write("];\n")

  nodes,values,triangles,extra_polys = info

  coords = num.array(nodes)

  x = coords[:,0]
  y = coords[:,1]

  writeMatlabVector("x", x)
  writeMatlabVector("y", y)
	  
  writeMatlabVector("node_data", num.array(values))

  tris = num.array([ [a+1,b+1,c+1] for (a,b,c) in triangles ])

  writeMatlabMatrix("tris", tris) 

  m_file.write("trisurf(tris, x, y, node_data)")





def _visualizationDriver(driver, filename, info):
  if driver == "gnuplot":
    _writeGnuplotFile(filename, info)
  elif driver == "matlab":
    _writeMatlabFile(filename, info)
  elif driver == "vtk":
    file_main, file_grid = filename
    _writeVtkFile(file_main, file_grid, info)
  else:
    raise RuntimeError, "invalid visualization driver: %s" % driver




def visualize(driver, filename, mesh, vector):
  _visualizationDriver(driver, filename, 
                       _compileInfo(mesh.elements(), vector))




def visualizeSeveralMeshes(driver, filename, offsets_meshes_and_vectors):
  nodes = []
  values = []
  triangles = []
  extra_polys = []
  for offset, mesh, vector in offsets_meshes_and_vectors:
    my_nodes,my_values,my_triangles,my_extra_polys = _compileInfo(mesh.elements(), vector)
    node_number_offset = len(nodes)

    def remap_node_numbers(poly_spec):
      return [[i+node_number_offset for i in node_number_list] 
              for node_number_list in poly_spec]

    triangles += remap_node_numbers(my_triangles)
    extra_polys += remap_node_numbers(my_extra_polys)
    values += my_values
    nodes += [node + offset for node in my_nodes]

  info = (nodes, values, triangles, extra_polys)
  _visualizationDriver(driver, filename, info)





def visualizePerElementData(name, visualizer, mesh, func_on_elements):
  raise RuntimeError, "broken"
  node_to_element_hash = {}
  for el in mesh.elements():
    for node in el.nodes():
      try:
        node_to_element_hash[node].append(el)
      except KeyError:
        node_to_element_hash[node] = [el]
  
  dof_count = mesh.dofManager().countDegreesOfFreedom()
  node_data = num.zeros((dof_count,), num.Float)
  for node_number in range(dof_count):
    node = mesh.dofManager().getDegreeOfFreedomIdentifier(node_number)
    node_data[node_number] = sum(map(func_on_elements, node_to_element_hash[node]))
  visualizer(name, mesh.dofManager(), mesh.elements(), node_data)




def writeGnuplotSparsityPattern(name, matrix):
  gnuplot_file = file(name, "w")
  for i,j in matrix.indices():
    gnuplot_file.write("%d %d\n" % (i,j))

