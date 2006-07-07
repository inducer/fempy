import pylinear.array as num
import mesh_function
import pytools




class VisualizationData:
    def __init__(self, start_nodes, start_values):
        self._Nodes = [node.Coordinates for node in start_nodes]
        self._Values = start_values
        self._Triangles = []
        self._ExtraPolygons = []
    
    def get_info_tuple(self):
        return self._Nodes, self._Values, self._Triangles, \
               self._ExtraPolygons

    def register_local_node(self, coordinates, data):
        idx = len(self._Nodes)
        self._Nodes.append(coordinates)
        self._Values.append(data)
        return idx

    def add_triangle(self, a, b, c):
        self._Triangles.append((a,b,c))
        
    def add_extra_polygon(self, nodes):
        self._ExtraPolygons.append(nodes)





def _write_gnuplot_file(name, info):
    gnuplot_file = file(name, "w")

    def write_node(node):
        gnuplot_file.write("%f %f %f\n" % (
            nodes[node][0],
            nodes[node][1],
            values[node]))
        
    nodes,values,triangles,extra_polys = info
    for tri_nodes in triangles:
        for node in tri_nodes:
            write_node(node)
            write_node(tri_nodes[0])
        gnuplot_file.write("\n\n")




def _write_vtk_file(name, name_grid, info):
    import pyvtk
    
    nodes,values,triangles,extra_polys = info
    
    # FIXME: temporary hack to work around VTK bug
    values2 = []
    for v in values:
        if abs(v) < 1e-38:
            values2.append(0)
        else:
            values2.append(v)

    my_nodes = [[node[0], node[1], 0.] for node in nodes]
      
    structure = pyvtk.PolyData(points = my_nodes, polygons = triangles)
    structure_grid = pyvtk.PolyData(points = my_nodes, polygons = extra_polys)
        
    pointdata = pyvtk. PointData(
        pyvtk. Scalars(values2, name="node_data", lookup_table = "default"))

    vtk = pyvtk.VtkData(structure, "FEM result triangles", pointdata)
    vtk.tofile(name, "ascii")

    vtk = pyvtk.VtkData(structure_grid, "FEM result grid", pointdata)
    vtk.tofile(name_grid, "ascii")



def _write_matlab_file(name, info):
    m_file = file(name, "w")

    def write_matlab_vector(name, data):
        (h,) = data.shape
        m_file.write("%s = [\n" % name)
        for i in range(0,h):
            if i != h-1:
                m_file.write("%f;\n" % data[ i ])
            else:
                m_file.write("%f" % data[ i ])
        m_file.write("];\n")

    def write_matlab_matrix(name, data):
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

    write_matlab_vector("x", x)
    write_matlab_vector("y", y)
	  
    write_matlab_vector("node_data", num.array(values))

    tris = num.array([ [a+1,b+1,c+1] for (a,b,c) in triangles ])
    
    write_matlab_matrix("tris", tris) 

    m_file.write("trisurf(tris, x, y, node_data)")





def _visualization_driver(driver, filename, info):
    if driver == "gnuplot":
        _write_gnuplot_file(filename, info)
    elif driver == "matlab":
        _write_matlab_file(filename, info)
    elif driver == "vtk":
        file_main, file_grid = filename
        _write_vtk_file(file_main, file_grid, info)
    else:
        raise RuntimeError, "invalid visualization driver: %s" % driver




def _get_visualization_data(mesh_func):
    dof_manager = mesh_func.mesh().dof_manager()

    nonu_lookup = pytools.reverse_dictionary(mesh_func.number_assignment())

    nodes = [nonu_lookup[nonu] for nonu in range(len(dof_manager))]
    node_values = list(mesh_func.vector())

    vis_data = VisualizationData(nodes, node_values)

    for el in mesh_func.mesh().elements():
        el.get_visualization_data(mesh_func, vis_data)
    return vis_data




def visualize(driver, filename, mesh_func):
    vis_data = _get_visualization_data(mesh_func)
    _visualization_driver(driver, filename, vis_data.get_info_tuple())




def visualize_several_meshes(driver, filename, offsets_and_mesh_functions):
    nodes = []
    values = []
    triangles = []
    extra_polys = []
    for offset, mesh_func in offsets_and_mesh_functions:
        vis_data = _get_visualization_data(mesh_func)
        my_nodes,my_values,my_triangles,my_extra_polys = vis_data.get_info_tuple()

        node_number_offset = len(nodes)

        def remap_node_numbers(poly_spec):
            return [[i+node_number_offset for i in node_number_list] 
                    for node_number_list in poly_spec]

        triangles += remap_node_numbers(my_triangles)
        extra_polys += remap_node_numbers(my_extra_polys)
        values += my_values
        nodes += [node + offset for node in my_nodes]

    info = (nodes, values, triangles, extra_polys)
    _visualization_driver(driver, filename, info)




def visualize_per_element_data(name, visualizer, mesh, func_on_elements):
    """Broken.
    """
    node_to_element_hash = {}
    for el in mesh.elements():
        for node in el.nodes():
            try:
                node_to_element_hash[node].append(el)
            except KeyError:
                node_to_element_hash[node] = [el]
  
    dof_count = mesh.dof_manager().countDegreesOfFreedom()
    node_data = num.zeros((dof_count,), num.Float)
    for node_number in range(dof_count):
        node = mesh.dof_manager().getDegreeOfFreedomIdentifier(node_number)
        node_data[node_number] = sum(map(func_on_elements, node_to_element_hash[node]))
    visualizer(name, mesh.dof_manager(), mesh.elements(), node_data)




def write_gnuplot_sparsity_pattern(name, matrix):
    gnuplot_file = file(name, "w")
    for i,j in matrix.indices():
        gnuplot_file.write("%d %d\n" % (i,j))


def write_gnuplot_mesh(mesh, filename):
    mesh_func = mesh_function.discretize_function(mesh, lambda x: 0)

    def write_node(node_number):
        gnuplot_file.write("%f %f\n" % (
            nodes[node_number][0],
            nodes[node_number][1] ))
        
    nodes,values,triangles,extra_polys = _get_visualization_data(mesh_func).get_info_tuple()

    gnuplot_file = file(filename, "w")
    for node_numbers in extra_polys:
        for node_number in node_numbers:
            write_node(node_number)
        write_node(node_numbers[0])
        gnuplot_file.write("\n\n")


def visualize_1d_mesh_function(mf, filename):
    gnuplot_file = file(filename, "w")
    nodes = [node for node in mf.mesh().dof_manager()]
    nodes.sort(lambda n1, n2: cmp(n1.Coordinates[0], n2.Coordinates[0]))
    for n in nodes:
        print >> gnuplot_file, "%f\t%f" % (n.Coordinates[0], mf[n])


def _visualize_2d_grid_data_v_t_k(major_grid, minor_grid, data, filename, 
                            scale_major = 1., scale_minor = 1.):
    import pyvtk

    points = [(scale_major * x,scale_minor * y,0) for x in major_grid for y in minor_grid]
    structure = pyvtk.StructuredGrid((len(minor_grid), len(major_grid), 1), points)
        
    flat_data = []
    for i in data:
        flat_data += i
    pointdata = pyvtk. PointData(
        pyvtk. Scalars(flat_data, name="node_data", lookup_table = "default"))

    vtk = pyvtk.VtkData(structure, "grid data", pointdata)
    vtk.tofile(filename, "ascii")


def _visualize_2d_grid_data_gnuplot(major_grid, minor_grid, data, filename, 
                                scale_major = 1., scale_minor = 1.):
    gpfile = file(filename, "w")

    def write_node(imajor, iminor):
        gpfile.write("%f\t%f\t%f\n" % (major_grid[imajor] * scale_major, 
                                         minor_grid[iminor] * scale_minor,
                                         data[imajor][iminor]))

    for ix in range(len(major_grid)-1):
        for iy in range(len(minor_grid)-1):
            write_node(ix, iy)
            write_node(ix+1, iy)
            write_node(ix+1, iy+1)
            write_node(ix, iy+1)
            gpfile.write("\n\n")
            



def visualize_2d_grid_data(driver, major_grid, minor_grid, data, filename, 
                        scale_major = 1., scale_minor = 1.):
    if driver == "vtk":
        _visualize_2d_grid_data_v_t_k(major_grid, minor_grid, data, filename,
                                scale_major, scale_minor)
    elif driver == "gnuplot":
        _visualize_2d_grid_data_gnuplot(major_grid, minor_grid, data, filename,
                                    scale_major, scale_minor)

