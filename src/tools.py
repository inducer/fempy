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




def indices_and(list):
  return zip(range(len(list)), list)




def distanceToLine(start_point, direction, point):
  # Ansatz: start_point + alpha * direction 
  # <start_point + alpha * direction - point, direction> = 0!
  alpha = - num.innerproduct(start_point - point, direction)/tools.norm2squared(direction)
  foot_point = start_point + alpha * direction
  return norm2(point - foot_point), alpha




def angleBetweenVectors(vec1, vec2):
  return math.acos(math.fabs(num.innerproduct(vec1, vec2) \
                             / norm2(vec1) \
                             / norm2(vec2)))




class tFunctionValueCache:
  def __init__(self, function):
    self.Function = function
    self.ResultMap = {}

  def __call__(self, arg):
    try:
      return self.ResultMap[arg]
    except KeyError:
      result = self.Function(arg)
      self.ResultMap[arg] = result
      return result




def sumOver(function, arguments):
  result = 0
  for i in arguments:
    result += function(i)
  return result




