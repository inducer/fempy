import pylinear.matrices as num
import stopwatch
import element
import mesh_function
import tools




class tVisualizationData:
    def __init__(self, start_nodes, start_values):
        self._Nodes = [node.Coordinates for node in start_nodes]
        self._Values = start_values
        self._Triangles = []
        self._ExtraPolygons = []
    
    def getInfoTuple(self):
        return self._Nodes, self._Values, self._Triangles, \
               self._ExtraPolygons

    def registerLocalNode(self, coordinates, data):
        idx = len(self._Nodes)
        self._Nodes.append(coordinates)
        self._Values.append(data)
        return idx

    def addTriangle(self, a, b, c):
        self._Triangles.append((a,b,c))
        
    def addExtraPolygon(self, nodes):
        self._ExtraPolygons.append(nodes)





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
    
    my_nodes = [[node[0], node[1], 0.] for node in nodes]
      
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




def _getVisualizationData(mesh_func):
    # also assign numbers for constrained nodes
    dof_manager = mesh_func.mesh().dofManager()
    new_number_assignment = element.assignNodeNumbers(dof_manager,
                                                      mesh_func.numberAssignment())

    my_mesh_func = mesh_func.copy(number_assignment = new_number_assignment)

    nonu_lookup = tools.reverseDictionary(new_number_assignment)

    nodes = [nonu_lookup[nonu] for nonu in range(len(dof_manager))]
    node_values = [mesh_func[node] for node in nodes]

    vis_data = tVisualizationData(nodes, node_values)

    for el in mesh_func.mesh().elements():
        el.getVisualizationData(my_mesh_func, vis_data)
    return vis_data




def visualize(driver, filename, mesh_func):
    vis_data = _getVisualizationData(mesh_func)
    _visualizationDriver(driver, filename, vis_data.getInfoTuple())




def visualizeSeveralMeshes(driver, filename, offsets_and_mesh_functions):
    nodes = []
    values = []
    triangles = []
    extra_polys = []
    for offset, mesh_func in offsets_and_mesh_functions:
        vis_data = _getVisualizationData(mesh_func)
        my_nodes,my_values,my_triangles,my_extra_polys = vis_data.getInfoTuple()

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
    """Broken.
    """
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

