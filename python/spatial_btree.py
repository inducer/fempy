import Numeric as num
import types




def doBoxesIntersect( (bl1,tr1), (bl2,tr2) ):
  (dimension,) = bl1.shape
  for i in range( 0, dimension ):
    if max( bl1[i], bl2[i] ) > min( tr1[i], tr2[i] ):
      return False
  return True




def getElementsBoundingBox( elements ):
  if len( elements ) == 0:
    raise RuntimeError, "Cannot get the bounding box of no elements."

  bboxes = [ el.boundingBox() for el in elements ]
  bottom_lefts = [ bl for bl,tr in bboxes ]
  top_rights = [ tr for bl,tr in bboxes ]
  return num.minimum.reduce( bottom_lefts ), num.maximum.reduce( top_rights )





def makeBuckets( bottom_left, top_right, allbuckets ):
  (dimensions,) = bottom_left.shape

  half = (top_right - bottom_left) / 2.
  def do( dimension, pos ):
    if dimension == dimensions:
      origin = bottom_left + pos * half 
      bucket = tSpatialBinaryTreeBucket( origin, origin + half )
      allbuckets.append( bucket )
      return bucket
    else:
      pos[ dimension ] = 0
      first = do( dimension + 1, pos )
      pos[ dimension ] = 1
      second = do( dimension + 1, pos )
      return [ first, second ]

  return do( 0, num.zeros( (dimensions,) ) )




class tSpatialBinaryTreeBucket:
  """This class represents one bucket in a spatial binary tree.
  It automatically decides whether it needs to create more subdivisions 
  beneath itself or not.
  """

  def __init__( self, bottom_left, top_right ):
    self.BottomLeft = bottom_left
    self.TopRight = top_right
    self.Center = (bottom_left + top_right) / 2
    
    # As long as Buckets is None, there are no subdivisions
    self.Buckets = None
    self.Elements = []
    
  def bottomLeft( self ):
    return self.BottomLeft

  def topRight( self ):
    return self.TopRight

  def addElement( self, element ):
    (dimensions,) = self.BottomLeft.shape
    if self.Buckets is None:
      # No subdivisions yet.
      if len( self.Elements ) > 8 * 2**dimensions:
	# Too many elements. Need to do subdivision.
	self.AllBuckets = []
	self.Buckets = makeBuckets( self.BottomLeft, self.TopRight, 
	  self.AllBuckets )

	for el in self.Elements:
	  self.insertIntoSubdivision( el )
	
	# Free up some memory
	del self.Elements

	self.insertIntoSubdivision( element )
      else:
	# Simple:
	self.Elements.append( element )
    else:
      self.insertIntoSubdivision( element )

  def insertIntoSubdivision( self, element ):
    # Surprisingly, this can't be done much faster.
    bbox = element.boundingBox()
    for bucket in self.AllBuckets:
      if doBoxesIntersect( ( bucket.bottomLeft(), bucket.topRight() ), bbox ):
	bucket.addElement( element )

  def findElement( self, point ):
    if self.Buckets:
      # We have subdivisions. Use them.
      (dimensions,) = point.shape
      bucket = self.Buckets
      for dim in range( 0, dimensions ):
	if point[ dim ] < self.Center[ dim ]:
	  bucket = bucket[ 0 ]
	else:
	  bucket = bucket[ 1 ]
      
      return bucket.findElement( point )
    else:
      # We don't. Perform linear search.
      for el in self.Elements:
	if el.isInElement( point ):
	  return el
      return None

  def visualize( self, file ):
    file.write( "%f %f\n" % ( self.BottomLeft[0], self.BottomLeft[1] ) );
    file.write( "%f %f\n" % ( self.TopRight[0], self.BottomLeft[1] ) );
    file.write( "%f %f\n" % ( self.TopRight[0], self.TopRight[1] ) );
    file.write( "%f %f\n" % ( self.BottomLeft[0], self.TopRight[1] ) );
    file.write( "%f %f\n\n" % ( self.BottomLeft[0], self.BottomLeft[1] ) );
    if self.Buckets:
      for i in self.AllBuckets:
	i.visualize( file )




def buildSpatialBinaryTree( elements ):
  bottom_left,top_right = getElementsBoundingBox( elements )
  bucket = tSpatialBinaryTreeBucket( bottom_left, top_right )
  for el in elements:
    bucket.addElement( el )
  return bucket




def buildElementFinder( elements ):
  bucket = buildSpatialBinaryTree( elements )
  return lambda point: bucket.findElement( point )
