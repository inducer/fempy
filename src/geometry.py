import math
import mesh
import expression_operators as eo




def getParallelogram(grid_vectors):
  """Returns a list of points (as tuples) that represent
  the parallelogram given by the parameters.
  """
  gv0 = grid_vectors[0] * 0.5
  gv1 = grid_vectors[1] * 0.5
  return [gv0+gv1, -gv0+gv1, -gv0-gv1, gv0-gv1]




def getCircle(radius, use_exact = True):
  sqrt2_inv = math.sqrt(2)/2 * radius
  r_squared = radius * radius

  return [
    # going counterclockwise.
    #
    #  T
    # L R ^
    #  B 

    # right
    mesh.tShapeGuide(0, [-sqrt2_inv, sqrt2_inv],
                (eo.POWER,(eo.MINUS,r_squared,(eo.POWER,(eo.VARIABLE,"t"),2)),0.5),
                use_exact_elements = use_exact),
    
    # top
    mesh.tShapeGuide(1, [sqrt2_inv, -sqrt2_inv],
                (eo.POWER,(eo.MINUS,r_squared,(eo.POWER,(eo.VARIABLE,"t"),2)),0.5),
                use_exact_elements = use_exact),

    # left
    mesh.tShapeGuide(0, [sqrt2_inv, -sqrt2_inv],
                (eo.NEG,(eo.POWER,(eo.MINUS,r_squared,(eo.POWER,(eo.VARIABLE,"t"),2)),0.5)),
                use_exact_elements = use_exact),

    # bottom
    mesh.tShapeGuide(1, [-sqrt2_inv,sqrt2_inv],
                (eo.NEG,(eo.POWER,(eo.MINUS,r_squared,(eo.POWER,(eo.VARIABLE,"t"),2)),0.5)),
                use_exact_elements = use_exact),
    ]




def getAnnulusGeometry(outer_radius, inner_radius, use_exact = True):
  return [mesh.tShapeSection(getCircle(outer_radius, use_exact), "dirichlet"),
          mesh.tShapeSection(getCircle(inner_radius, use_exact), "dirichlet")]
