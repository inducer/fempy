import math
import pylinear.matrices as num
import pylinear.linear_algebra as la
import pylinear.matrix_tools as mtools




def estimateOrderOfConvergence(abscissae, errors):
  """Assuming that abscissae and errors are connected by a law of the form

  error = constant * abscissa ^ (-order),

  this function finds, in a least-squares sense, the best approximation of
  constant and order for the given data set. It returns a tuple (constant, order).
  Both inputs must be PyLinear vectors.
  """

  assert len(abscissae) == len(errors)
  coefficients = mtools.fit_polynomial(num.log10(abscissae), num.log10(errors), 1)
  return 10**coefficients[0], -coefficients[1]


  




class tEOCRecorder:
  def __init__(self):
    self.History = []

  def addDataPoint(self, abscissa, error):
    self.History.append((abscissa, error))

  def estimateOrderOfConvergence(self, gliding_mean = None):
    abscissae = num.array([ a for a,e in self.History ])
    errors = num.array([ e for a,e in self.History ])

    size = len(abscissae)
    if gliding_mean is None:
      gliding_mean = size

    data_points = size - gliding_mean + 1
    result = num.zeros((data_points, 2), num.Float)
    for i in range(data_points):
      result[i,0], result[i,1] = estimateOrderOfConvergence(
        abscissae[i:i+gliding_mean], errors[i:i+gliding_mean])
    return result

  def writeGnuplotFile(self, filename):
    outfile = file(filename, "w")
    for absc, err in self.History:
      outfile.write("%f %f\n" % (math.log10(absc), math.log10(err)))
    result = self.estimateOrderOfConvergence()
    const = result[0,0]
    order = result[0,1]
    outfile.write("\n")
    for absc, err in self.History:
      outfile.write("%f %f\n" % (math.log10(absc), math.log10(const * absc**(-order))))

