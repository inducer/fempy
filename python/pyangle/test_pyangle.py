import pyangle
import math

segments = 50

points = [ (1,0),(1,1),(-1,1),(-1,-1),(1,-1),(1,0) ]

for i in range( 0, segments + 1 ):
  angle = i * 2 * math.pi / segments
  points.append( ( 0.5 * math.cos( angle ), 0.5 * math.sin( angle ) ) )

def needsRefinement( vert_origin, vert_destination, vert_apex, area ):
  bary_x = ( vert_origin.x() + vert_destination.x() + vert_apex.x() ) / 3
  bary_y = ( vert_origin.y() + vert_destination.y() + vert_apex.y() ) / 3

  dist_center = math.sqrt( bary_x**2 + bary_y**2 )
  max_area = math.fabs( 0.002 * (dist_center-0.5) ) + 0.0001
  return area > max_area

out_p = pyangle.triangulateArea( points, [ (0,0) ], verbose = True, 
    refinement_func = needsRefinement )

pyangle.writeGnuplotMesh( "+tris.dat", out_p )
