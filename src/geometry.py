import math
import mesh
import pymbolic




def get_parallelogram(grid_vectors):
    """Returns a list of points (as tuples) that represent
    the parallelogram given by the parameters.
    """
    gv0 = grid_vectors[0] * 0.5
    gv1 = grid_vectors[1] * 0.5
    return [gv0+gv1, -gv0+gv1, -gv0-gv1, gv0-gv1]




def get_circle(radius, use_exact=True):
    sqrt2_inv = math.sqrt(2)/2 * radius
    r_squared = pymbolic.const(radius * radius)
    v_t = pymbolic.var("t")

    return [
        # going counterclockwise.
        #
        #  T
        # L R ^
        #  B 

        # right
        mesh.ShapeGuide(0, [-sqrt2_inv, sqrt2_inv],
                        (r_squared - v_t**2)**0.5,
                        use_exact_elements = use_exact),
    
        # top
        mesh.ShapeGuide(1, [sqrt2_inv, -sqrt2_inv],
                        (r_squared - v_t**2)**0.5,
                        use_exact_elements = use_exact),

        # left
        mesh.ShapeGuide(0, [sqrt2_inv, -sqrt2_inv],
                        -(r_squared - v_t**2)**0.5,
                        use_exact_elements = use_exact),
        
        # bottom
        mesh.ShapeGuide(1, [-sqrt2_inv,sqrt2_inv],
                        -(r_squared - v_t**2)**0.5,
                        use_exact_elements = use_exact),
        ]




def get_annulus_geometry(outer_radius, inner_radius, use_exact=True):
    return [mesh.ShapeSection(get_circle(outer_radius, use_exact), "dirichlet"),
            mesh.ShapeSection(get_circle(inner_radius, use_exact), "dirichlet")]
