import pylinear.matrices as num
import math

# Formulae from Schwarz, p. 121f.




def _doIntegration(locations, weights, f):
  result = 0
  for weight, location in zip(weights, locations):
    result += weight * f(location)
  return result




INTEGRATION_TRI_WEIGHTS = num.array([
    0.1125, # barycenter

    (155.+math.sqrt(15))/2400.,
    (155.+math.sqrt(15))/2400.,
    (155.+math.sqrt(15))/2400.,

    (155.-math.sqrt(15))/2400.,
    (155.-math.sqrt(15))/2400.,
    (155.-math.sqrt(15))/2400. ])

INTEGRATION_TRI_LOCATIONS = [
    num.array([ 1/3.,1/3., ]),

    num.array([ (6.+math.sqrt(15))/21., (6.+math.sqrt(15))/21. ]),
    num.array([ (9.-2.*math.sqrt(15))/21., (6.+math.sqrt(15))/21. ]),
    num.array([ (6.+math.sqrt(15))/21., (9.-2.*math.sqrt(15))/21. ]),

    num.array([ (6.-math.sqrt(15))/21., (6.-math.sqrt(15))/21. ]),
    num.array([ (9.+2.*math.sqrt(15))/21., (6.-math.sqrt(15))/21. ]),
    num.array([ (6.-math.sqrt(15))/21., (9.+2.*math.sqrt(15))/21. ]) ]

def integrateOnUnitTriangle(f):
  """This function returns the value of 

  \int_T f(x,y) d(x,y)

  where T is the right triangle with the cornes (0,0),(1,0),(0,1)."""
  return _doIntegration(INTEGRATION_TRI_LOCATIONS, INTEGRATION_TRI_WEIGHTS, f)

def integrateOnUnitTriangleWithOrder1(f):
  return 0.5 * f(num.array([1/3., 1/3.]))




  
INTEGRATION_LINE_LOCATIONS_3 = num.array([
    0.1127016654,
    0.5,
    0.8872983346 ])

INTEGRATION_LINE_WEIGHTS_3 = num.array([
    5./18.,
    8./18.,
    5./18. ])

def integrateOnUnitInterval3(f):
  """This function returns the value of

  \int_0^1 f(x) dx.
  
  It uses a Gaussian quadrature of order 3."""
  return _doIntegration(INTEGRATION_LINE_LOCATIONS_3, INTEGRATION_LINE_WEIGHTS_3, f)




INTEGRATION_LINE_LOCATIONS_4 = num.array([
    0.0694318442,
    0.3300094782,
    0.6699905218,
    0.9305681558 ])

INTEGRATION_LINE_WEIGHTS_4 = num.array([
    0.1739274226,
    0.3260725774,
    0.3260725774,
    0.1739274226 ])

def integrateOnUnitInterval4(f):
  """This function returns the value of

  \int_0^1 f(x) dx.
  
  It uses a Gaussian quadrature of order 4."""
  return _doIntegration(INTEGRATION_LINE_LOCATIONS_4, INTEGRATION_LINE_WEIGHTS_4, f)
