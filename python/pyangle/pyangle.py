import pyanglemetal
from pyanglemetal import tTriangulationParameters




def triangulate( in_parameters, verbose = False,
    refinement_func = None):
  opts = "pzq"
  if verbose:
    opts += "V"
  else:
    opts += "Q"
  if refinement_func is not None:
    opts += "u"
  out_parameters = tTriangulationParameters()
  pyanglemetal.triangulate( 
      opts, in_parameters, out_parameters, 
      tTriangulationParameters(),
      refinement_func )
  return out_parameters




def triangulateArea( points, hole_starts = [], verbose = False, refinement_func = None ):
  input_p = tTriangulationParameters()
  pointcount = len( points )

  input_p.Points.setSize( pointcount )
  input_p.Segments.setSize( pointcount )

  i = 0
  for pt in points:
    input_p.Points.setSub( i, 0, pt[0] )
    input_p.Points.setSub( i, 1, pt[1] )
    i += 1

  for i in range( 0, pointcount ):
    input_p.Segments.setSub( i, 0, i )
    input_p.Segments.setSub( i, 1, (i+1) % pointcount )

  holecount = len( hole_starts )
  input_p.Holes.setSize( len( hole_starts ) )
  for i in range( 0, holecount ):
    input_p.Holes.setSub( i, 0, hole_starts[i][0] )
    input_p.Holes.setSub( i, 1, hole_starts[i][1] )

  return triangulate( input_p, verbose, refinement_func )




def writeGnuplotMesh( gnuplot_file, out_p ):
  gp_file = file( gnuplot_file, "w" )

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

