import math
import fempy.mesh as mesh




def getParallelogram(edge_length = 1, x_skew = 0, y_skew = 0):
  """Returns a list of points (as tuples) that represent
  the parallelogram given by the parameters.
  The order of the points is as follows:

    10
    23
  """
  a = edge_length / 2.
  xs = x_skew / 2.
  ys = y_skew / 2.
  return [(a-xs,-a+ys), (a+xs,a+ys), (-a+xs,a-ys), (-a-xs,-a-ys)]

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
                ("**",("-",r_squared,("**",("variable","t"),2)),0.5),
                use_exact_elements = use_exact),
    
    # top
    mesh.tShapeGuide(1, [sqrt2_inv, -sqrt2_inv],
                ("**",("-",r_squared,("**",("variable","t"),2)),0.5),
                use_exact_elements = use_exact),

    # left
    mesh.tShapeGuide(0, [sqrt2_inv, -sqrt2_inv],
                ("-",("**",("-",r_squared,("**",("variable","t"),2)),0.5)),
                use_exact_elements = use_exact),

    # bottom
    mesh.tShapeGuide(1, [-sqrt2_inv,sqrt2_inv],
                ("-",("**",("-",r_squared,("**",("variable","t"),2)),0.5)),
                use_exact_elements = use_exact),
    ]

def getUnitCellGeometry(edge_length, segments = 50, inner_factor = 0.3, use_exact = True):
  return [mesh.tShapeSection(getParallelogram(edge_length), True),
          mesh.tShapeSection(getCircle(edge_length * inner_factor, use_exact), False)]

def getAnnulusGeometry(outer_radius, inner_radius, use_exact = True):
  return [mesh.tShapeSection(getCircle(outer_radius, use_exact), True),
          mesh.tShapeSection(getCircle(inner_radius, use_exact), True)]
