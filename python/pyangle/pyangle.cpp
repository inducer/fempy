#include "triangle.h"
#include <boost/python.hpp>
#include <vector>
#include <stdexcept>
#include <iostream>




using namespace boost::python;
using namespace std;




class tSizeChangeNotifier;




class tSizeChangeNotificationReceiver
{
  public:
    virtual void notifySizeChange( tSizeChangeNotifier *master, unsigned size ) = 0;
};




class tSizeChangeNotifier
{
    typedef vector<tSizeChangeNotificationReceiver *> tNotificationReceiverList;
    tNotificationReceiverList NotificationReceivers;

  public:
    virtual unsigned size() = 0;
    virtual void setSize( unsigned size )
    {
      tNotificationReceiverList::iterator first = NotificationReceivers.begin(),
      last = NotificationReceivers.end();
      while ( first != last )
	(*first++)->notifySizeChange( this, size );
    }

    void registerForNotification( tSizeChangeNotificationReceiver *rec )
    {
      NotificationReceivers.push_back( rec );
    }

    void unregisterForNotification( tSizeChangeNotificationReceiver *rec )
    {
      tNotificationReceiverList::iterator first = NotificationReceivers.begin(),
      last = NotificationReceivers.end();
      while ( first != last )
      {
	if ( rec == *first )
	{
	  NotificationReceivers.erase( first );
	  return;
	}
	first++;
      }
    }
};




class tVertex : public boost::noncopyable
{
  public:
    REAL	*Data;

  public:
    tVertex( REAL *data )
    : Data( data )
    {
    }

    REAL x() { return Data[ 0 ]; }
    REAL y() { return Data[ 1 ]; }
};




template<class ElementT> 
class tForeignArray : public tSizeChangeNotifier, public tSizeChangeNotificationReceiver,
  public boost::noncopyable
{
    ElementT	*&Contents;
    int		&NumberOf;
    unsigned	Unit;
    tSizeChangeNotifier *SlaveTo;

  public:
    tForeignArray( ElementT *&cts, int &number_of, unsigned unit = 1, tSizeChangeNotifier *slave_to = NULL )
      : Contents( cts ), NumberOf( number_of ), Unit( unit ), SlaveTo( slave_to )
    {
      Contents = NULL;
      if ( SlaveTo )
      {
	SlaveTo->registerForNotification( this );
	setSizeInternal( SlaveTo->size() );
      }
      else
	setSize( 0 );
    }

    ~tForeignArray()
    {
      if ( SlaveTo )
	SlaveTo->unregisterForNotification( this );

      if ( Contents )
	free( Contents );
      Contents = NULL;

      if ( !SlaveTo )
	NumberOf = 0;
    }

    unsigned size()
    {
      return NumberOf;
    }

    void setSize( unsigned size )
    {
      if ( SlaveTo )
	throw runtime_error( "sizes of slave arrays cannot be changed" );
      else
	setSizeInternal( size );
    }

    void notifySizeChange( tSizeChangeNotifier *master, unsigned size )
    {
      setSizeInternal( size );
    }

    void setSizeInternal( unsigned size )
    {
      if ( !SlaveTo )
	NumberOf = size;
      
      if ( Contents != NULL )
	free( Contents );

      if ( size == 0 || Unit == 0 )
	Contents = NULL;
      else
      {
	Contents = (ElementT *) malloc( sizeof( ElementT ) * Unit * size );
	if ( Contents == NULL )
	  throw bad_alloc();
      }

      tSizeChangeNotifier::setSize( size );
    }

    void setUnit( unsigned unit )
    {
      Unit = unit;
      setSize( NumberOf );
    }

    void set( unsigned index, ElementT value )
    {
      if ( index >= NumberOf * Unit )
	throw runtime_error( "index out of bounds" );
      Contents[ index ] = value;
    }

    void setSub( unsigned index, unsigned sub_index, ElementT value )
    {
      set( index * Unit + sub_index, value );
    }

    ElementT get( unsigned index )
    {
      if ( index >= NumberOf * Unit )
	throw runtime_error( "index out of bounds" );
      return Contents[ index ];
    }

    ElementT getSub( unsigned index, unsigned sub_index )
    {
      return get( index * Unit + sub_index );
    }
};




struct tTriangulationParameters : public triangulateio, public boost::noncopyable
{
  public:
    tForeignArray<REAL>		Points; // in/out
    tForeignArray<REAL>		PointAttributes; // in/out
    tForeignArray<int>		PointMarkers; // in/out

    tForeignArray<int>		Triangles; // in/out
    tForeignArray<REAL>		TriangleAttributes; // in/out
    tForeignArray<REAL>		TriangleAreas; // in only
    tForeignArray<int>		Neighbors; // out only

    tForeignArray<int>		Segments; // in/out
    tForeignArray<int>		SegmentMarkers; // in/out
    
    tForeignArray<REAL>		Holes; // in only

    tForeignArray<REAL>		Regions; // in only

    tForeignArray<int>		Edges; // out only
    tForeignArray<int>		EdgeMarkers; // out only
    tForeignArray<REAL>		Normals; // out only

  public:
    tTriangulationParameters()
      : Points( pointlist, numberofpoints, 2 ),
        PointAttributes( pointattributelist, numberofpoints, 0, &Points ),
	PointMarkers( pointmarkerlist, numberofpoints, 1, &Points ),

	Triangles( trianglelist, numberoftriangles, 3 ),
	TriangleAttributes( triangleattributelist, numberoftriangles, 0, &Triangles ),
	TriangleAreas( trianglearealist, numberoftriangles, 1, &Triangles ),
	Neighbors( neighborlist, numberoftriangles, 3, &Triangles ),

	Segments( segmentlist, numberofsegments, 2 ),
	SegmentMarkers( segmentmarkerlist, numberofsegments, 1, &Segments ),

	Holes( holelist, numberofholes, 2 ),

	Regions( regionlist, numberofregions, 4 ),

	Edges( edgelist, numberofedges, 2 ),
	EdgeMarkers( edgemarkerlist, numberofedges, 1, &Edges ),
	Normals( normlist, numberofedges, 2, &Edges )
    {
      numberofpointattributes = 0;
      numberofcorners = 3;
      numberoftriangleattributes = 0;
    }

    void setNumberOfPointAttributes( unsigned attrs )
    {
      PointAttributes.setUnit( attrs );
      numberofpointattributes = attrs;
    }

    void setNumberOfTriangleAttributes( unsigned attrs )
    {
      TriangleAttributes.setUnit( attrs );
      numberoftriangleattributes = attrs;
    }
};




object RefinementFunction;




int triunsuitable( vertex triorg, vertex tridest, vertex triapex, REAL area )
{
  // return 1 if triangle is too large, 0 otherwise
  try
  {
    tVertex org( triorg );
    tVertex dest( tridest );
    tVertex apex( triapex );
    return extract<bool>( RefinementFunction( 
	  boost::ref( org ), boost::ref( dest ), boost::ref( apex ), area ) );
  }
  catch ( exception &ex )
  {
    cerr 
      << "*** Oops. Your Python refinement function raised an exception." << endl
      << "*** " << ex.what() << endl
      << "*** Sorry, we can't continue." << endl;
    abort();
  }
  catch (...)
  {
    cerr 
      << "*** Oops. Your Python refinement function raised an exception." << endl
      << "*** Sorry, we can't continue." << endl;
    abort();
  }
}




void triangulateWrapper(char *options, tTriangulationParameters &in, 
    tTriangulationParameters &out,
    tTriangulationParameters &voronoi,
    object refinement_func)
{
  RefinementFunction = refinement_func;

  REAL u[] = {
    1,2
  };
  tVertex vert( u );
  RefinementFunction( boost::ref( vert ), boost::ref( vert ), boost::ref( vert ), 5 );
  
  triangulate( options, &in, &out, &voronoi );

  RefinementFunction = object(); // i.e. None
  out.holelist = NULL;
  out.numberofholes = 0;

  out.regionlist = NULL;
  out.numberofregions = 0;
}




BOOST_PYTHON_MODULE(pyanglemetal)
{
  def( "triangulate", triangulateWrapper );
  class_<tTriangulationParameters, bases<>, tTriangulationParameters, boost::noncopyable>
    ( "tTriangulationParameters" )
    .def_readonly( "Points", &tTriangulationParameters::Points )
    .def_readonly( "PointAttributes", &tTriangulationParameters::Points )
    .def_readonly( "PointMarkers", &tTriangulationParameters::PointMarkers )

    .def_readonly( "Triangles", &tTriangulationParameters::Triangles )
    .def_readonly( "TriangleAttributes", &tTriangulationParameters::TriangleAttributes )
    .def_readonly( "TriangleAreas", &tTriangulationParameters::TriangleAreas )
    .def_readonly( "Neighbors", &tTriangulationParameters::Neighbors )

    .def_readonly( "Segments", &tTriangulationParameters::Segments )
    .def_readonly( "SegmentMarkers", &tTriangulationParameters::SegmentMarkers )

    .def_readonly( "Holes", &tTriangulationParameters::Holes )

    .def_readonly( "Regions", &tTriangulationParameters::Regions )

    .def_readonly( "Edges", &tTriangulationParameters::Edges )
    .def_readonly( "EdgeMarkers", &tTriangulationParameters::EdgeMarkers )

    .def_readonly( "Normals", &tTriangulationParameters::Normals );
  
  class_<tForeignArray<REAL>, bases<>, tForeignArray<REAL>, boost::noncopyable >
    ( "tRealArray", no_init )
    .def( "size", &tForeignArray<REAL>::size )
    .def( "setSize", &tForeignArray<REAL>::setSize )
    .def( "set", &tForeignArray<REAL>::set )
    .def( "setSub", &tForeignArray<REAL>::setSub )
    .def( "get", &tForeignArray<REAL>::get )
    .def( "getSub", &tForeignArray<REAL>::getSub );
	
  class_<tForeignArray<int>, bases<>, tForeignArray<int>, boost::noncopyable >
    ( "tIntArray", no_init )
    .def( "size", &tForeignArray<int>::size )
    .def( "setSize", &tForeignArray<int>::setSize )
    .def( "set", &tForeignArray<int>::set )
    .def( "setSub", &tForeignArray<int>::setSub )
    .def( "get", &tForeignArray<int>::get )
    .def( "getSub", &tForeignArray<int>::getSub );

  class_<tVertex, bases<>, tVertex, boost::noncopyable>( "tVertex", no_init )
    .def( "x", &tVertex::x )
    .def( "y", &tVertex::y );
}
