import Numeric as num




def writeGnuplotFile( name, dof_manager, elements, solution ):
  gnuplot_file = file( name, "w" )

  def writeNode( node ):
    gnuplot_file.write( "%f %f %f\n" % (
	node.coordinates()[0],
	node.coordinates()[1],
	solution[ dof_manager.getDegreeOfFreedomNumber( node ) ] ) )

  for element in elements:
    for node in element.nodes(): writeNode( node )
    writeNode( element.nodes()[0] )
    gnuplot_file.write( "\n\n" )




def writeVtkFile( name, dof_manager, elements, solution ):
  import pyvtk

  dof_count = dof_manager.countDegreesOfFreedom()

  points_list = []
  for dof in range( 0, dof_count ):
    points_list.append( vector2tuple (
	dof_manager.getDegreeOfFreedomIdentifier( dof ).coordinates() ) )

  polygon_list = []
  for el in elements:
    polygon_list.append( el.nodeNumbers() )

  structure = pyvtk.PolyData( points=points_list, polygons=polygon_list)

  solution_list = []
  for i in solution:
    solution_list.append( i )

  pointdata = pyvtk. PointData(
      pyvtk. Scalars(solution_list, name="solution", lookup_table = "default") )

  vtk = pyvtk.VtkData( structure, "FEM result", pointdata )
  vtk.tofile( name, "ascii" )




def writeMatlabFile( name, dof_manager, elements, solution ):
  m_file = file( name, "w" )

  def writeMatlabVector( name, data ):
    (h,) = data.shape
    m_file.write( "%s = [\n" % name )
    for i in range( 0,h ):
      if i != h-1:
	m_file.write( "%f;\n" % data[ i ] )
      else:
	m_file.write( "%f" % data[ i ] )
    m_file.write( "];\n" )

  def writeMatlabMatrix( name, data ):
    (h,w) = data.shape
    m_file.write( "%s = [\n" % name )
    for i in range( 0,h ):
      for j in range( 0,w ):
	if j != w-1:
	  m_file.write( "%f," % data[ i,j ] )
	else:
	  if i != h-1:
	    m_file.write( "%f;\n" % data[ i,j ] )
	  else:
	    m_file.write( "%f" % data[ i,j ] )
    m_file.write( "];\n" )

  writeMatlabVector( "solution", solution )

  dof_count = dof_manager.countDegreesOfFreedom()
  
  x = num.zeros( (dof_count,), num.Float )
  y = num.zeros( (dof_count,), num.Float )

  for i in range( 0, dof_count ):
    coords = dof_manager.getDegreeOfFreedomIdentifier( i ).coordinates()
    x[ i ] = coords[ 0 ]
    y[ i ] = coords[ 1 ]

  writeMatlabVector( "x", x )
  writeMatlabVector( "y", y )
	  
  tris = num.zeros( (len(elements),3), num.Float )

  el_index = 0
  for el in elements:
    tris[ el_index ] = num.array( el.nodeNumbers() )
    el_index += 1

  tris += num.ones( tris.shape, num.Float )

  writeMatlabMatrix( "tris", tris ) 

  m_file.write( "trisurf( tris, x, y, solution )" )





