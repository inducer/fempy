import pylinear.matrices as num
import math




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
    f_value = x_start
    x_start -= f_value / fprime(x_start)
    if math.fabs(f_value) < tolerance:
      return x_start
  raise RuntimeError, "Newton iteration failed, a zero was not found"





