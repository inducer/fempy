import pylinear.matrices as num
import pylinear.linear_algebra as la
import math




class tReference:
  def __init__( self, value ):
    self.V = value
  def get( self ):
    return self.V
  def set( self, value ):
    self.V = value




def flatten(list):
  result = []
  for i in list:
    result.extend(i)
  return result




def writeMatrixAsCSV(filename, matrix):
  mat_file = file(filename, "w")
  h,w = matrix.shape
  for row in range(0, h):
    for column in range(0, w):
      mat_file.write("%f," % matrix[ row, column ])
    mat_file.write("\n")




def norm2squared(vector):
  return num.innerproduct(vector, vector)




def norm2(vector):
  return math.sqrt(num.innerproduct(vector, vector))




def sequence2list(sequence):
  result = []
  for i in sequence:
    result.append(i)
  return result




def vector2tuple(vector):
  size, = vector.shape
  if size == 2:
    return vector[ 0 ], vector[ 1 ], 0
  else:
    return vector[ 0 ], vector[ 1 ], vector[ 2 ]




def findZeroByNewton(f, fprime, x_start, tolerance = 1e-12, maxit = 10):
  it = 0
  while it < maxit:
    it += 1
    f_value = f(x_start)
    x_start -= f_value / fprime(x_start)
    if math.fabs(f_value) < tolerance:
      return x_start
  raise RuntimeError, "Newton iteration failed, a zero was not found"




def findVectorZeroByNewton(f, fprime, x_start, tolerance = 1e-12, maxit = 10):
  it = 0
  while it < maxit:
    it += 1
    f_value = f(x_start)
    x_start -= num.matrixmultiply(la.inverse(fprime(x_start)), f_value)
    if norm2(f_value) < tolerance:
      return x_start
  raise RuntimeError, "Newton iteration failed, a zero was not found"




def average(sequence):
  return sum(sequence)/float(len(sequence))



def allEqual(sequence):
  item = sequence[0]
  for i in sequence[1:]:
    if i != item:
      return False
  return True




def decorate(function, list):
  return map(lambda x: (x, function(x)), list)




def partition(criterion, list):
  part_true = []
  part_false = []
  for i in list:
    if criterion(i):
      part_true.append(i)
    else:
      part_false.append(i)
  return part_true, part_false
