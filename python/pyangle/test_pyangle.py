import pyangle
import math

in_p = pyangle.tTriangulationParameters()
out_p = pyangle.tTriangulationParameters()

segments = 400
in_p.Segments.setSize( segments )
in_p.Points.setSize( segments )

segs = in_p.Segments
pts = in_p.Points

for i in range( 0, segments ):
  angle = i * 2 * math.pi / segments
  pts.setSub( i, 0,  math.cos( angle ) )
  pts.setSub( i, 1, math.sin( angle ) )

  segs.setSub( i, 0, i )
  segs.setSub( i, 1, (i + 1) % segments )

def needsRefinement( vert_origin, vert_destination, vert_apex, area ):
  bary_x = ( vert_origin.x() + vert_destination.x() + vert_apex.x() ) / 3
  bary_y = ( vert_origin.y() + vert_destination.y() + vert_apex.y() ) / 3

  dist_center = math.sqrt( bary_x**2 + bary_y**2 )
  max_area = (dist_center-0.5)**4 + 0.0001
  return area > max_area

pyangle.triangulate( in_p, out_p, verbose = True, 
    refinement_func = needsRefinement )

gp_file = file( "+tris.dat", "w" )

pts = out_p.Points
tris = out_p.Triangles

for tri in range( out_p.Triangles.size() ):
  for pt in range( 0, 3 ):
    gp_file.write( "%f %f\n" % (
      pts.getSub( tris.getSub( tri, pt ), 0 ),
      pts.getSub( tris.getSub( tri, pt ), 1 ) ) )
  gp_file.write( "%f %f\n" % (
    pts.getSub( tris.getSub( tri, 0 ), 0 ),
    pts.getSub( tris.getSub( tri, 0 ), 1 ) ) )
  gp_file.write( "\n" )
