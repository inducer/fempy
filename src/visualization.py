import pylinear.matrices as num




class tVisualizationData:
  """This object represents the visualization for one finite element.
  by specifying a number of triangles and the values of the solution
  function at their nodes. This is what the data mean:

  Nodes is a list of tNodes.

  NodeValues specifies the values of the solution function at the
  nodes given by NodeNumbers, in the same order.

  Triangles is a list of 3-tuples of node numbers. To refer to the
  above NodeNumbers array, use regular-zero based indexing into that
  array. But an element may introduce extra nodes into its 
  visualization that are not necessarily mentioned in the DOF manager.
  These are given in the Extra* arrays and referred to by negative
  indices in Triangles. (where -1 -> 0, -2 -> 1 and so on)

  ExtraNodeCoordinates is an array of Numeric arrays of (n,) shape
  (i.e. vectors), giving the coordinates of the extra nodes.

  ExtraNodeValues gives the values of the solution function at the
  extra nodes.
  """

  def __init__(self, nodes, node_values, triangles,
      extra_node_coordinates = [], extra_node_values = []):
    self.Nodes = nodes
    self.NodeValues = node_values
    self.Triangles = triangles
    self.ExtraNodeCoordinates = extra_node_coordinates
    self.ExtraNodeValues = extra_node_values




def compileInfo(dof_manager, elements, solution):
  node_number_map = {}
  nodes = []
  values = []
  tris = []
  for el in elements:
    data = el.visualizationData(solution)
    local_nodes = {}

    index = 0
    for node in data.Nodes:
      node_number = dof_manager.getDegreeOfFreedomNumber(node)
      if node_number in node_number_map:
	local_nodes[ index ] = node_number_map[ node_number ]
      else:
	local_nodes[ index ] = len(nodes)
	node_number_map[ node ] = len(nodes)
	nodes.append(node.coordinates())
	values.append(data.NodeValues[ index ])
      index += 1

    index = -1
    for coords in data.ExtraNodeCoordinates:
      local_nodes[ index ] = len(nodes)
      nodes.append(coords)
      values.append(data.ExtraNodeValues[ -(index+1) ])
      index -= 1

    for na,nb,nc in data.Triangles:
      tris.append((local_nodes[na], local_nodes[nb], local_nodes[nc]))

  return nodes, values, tris



def writeGnuplotFile(name, dof_manager, elements, solution):
  gnuplot_file = file(name, "w")

  def writeNode(node):
    gnuplot_file.write("%f %f %f\n" % (
	nodes[node][0],
	nodes[node][1],
	values[node]))

  nodes,values,triangles = compileInfo(dof_manager, elements, solution)
  for tri_nodes in triangles:
    for node in tri_nodes:
      writeNode(node)
    writeNode(tri_nodes[0])
    gnuplot_file.write("\n\n")




def writeVtkFile(name, dof_manager, elements, solution):
  import pyvtk

  nodes,values,triangles = compileInfo(dof_manager, elements, solution)

  my_nodes = []
  for node in nodes:
    my_nodes.append([ node[0], node[1], 0 ])

  structure = pyvtk.PolyData(points = my_nodes, polygons = triangles)

  pointdata = pyvtk. PointData(
      pyvtk. Scalars(values, name="solution", lookup_table = "default"))

  vtk = pyvtk.VtkData(structure, "FEM result", pointdata)
  vtk.tofile(name, "ascii")




def writeMatlabFile(name, dof_manager, elements, solution):
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

  nodes,values,triangles = compileInfo(dof_manager, elements, solution)

  coords = num.array(nodes)

  x = coords[:,0]
  y = coords[:,1]

  writeMatlabVector("x", x)
  writeMatlabVector("y", y)
	  
  writeMatlabVector("solution", num.array(values))

  tris = num.array([ [a+1,b+1,c+1] for (a,b,c) in triangles ])

  writeMatlabMatrix("tris", tris) 

  m_file.write("trisurf(tris, x, y, solution)")





def writeGnuplotSparsityPattern(name, matrix):
  gnuplot_file = file(name, "w")
  for i,j in matrix.indices():
    gnuplot_file.write("%d %d\n" % (i,j))

