import Numeric as num
import math




def flatten( list ):
  result = []
  for i in list:
    result.extend( i )
  return result




def writeSymmetricMatrixAsCSV( filename, matrix ):
  mat_file = file( filename, "w" )
  h,w = matrix.shape
  for row in range( 0, h ):
    for column in range( 0, row ):
      mat_file.write( "%f," % matrix[ row, column ] )
    for column in range( row, w ):
      mat_file.write( "%f," % matrix[ column, row ] )
    mat_file.write( "\n" )




def norm2( vector ):
  return math.sqrt( num.dot( vector, vector ) )




def sequence2list( sequence ):
  result = []
  for i in sequence:
    result.append( i )
  return result




def vector2tuple( vector ):
  size, = vector.shape
  if size == 2:
    return vector[ 0 ], vector[ 1 ], 0
  else:
    return vector[ 0 ], vector[ 1 ], vector[ 2 ]




