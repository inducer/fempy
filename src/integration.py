import Numeric as num
import math

# Formulae from Schwarz, p. 121f.




def doIntegration( locations, weights, f ):
  h,w = locations.shape
  f_values = num.zeros( (h,), num.Float )
  for i in range(0,h):
    f_values[ i ] = f( locations[ i ] )

  return num.dot( weights, f_values )




INTEGRATION_TRI_WEIGHTS = num.array( [
    0.1125, # barycenter

    (155.+math.sqrt(15))/2400.,
    (155.+math.sqrt(15))/2400.,
    (155.+math.sqrt(15))/2400.,

    (155.-math.sqrt(15))/2400.,
    (155.-math.sqrt(15))/2400.,
    (155.-math.sqrt(15))/2400. ] )

INTEGRATION_TRI_LOCATIONS = num.array( [
    [ 1/3.,1/3., ],

    [ (6.+math.sqrt(15))/21., (6.+math.sqrt(15))/21. ],
    [ (9.-2.*math.sqrt(15))/21., (6.+math.sqrt(15))/21. ],
    [ (6.+math.sqrt(15))/21., (9.-2.*math.sqrt(15))/21. ],

    [ (6.-math.sqrt(15))/21., (6.-math.sqrt(15))/21. ],
    [ (9.+2.*math.sqrt(15))/21., (6.-math.sqrt(15))/21. ],
    [ (6.-math.sqrt(15))/21., (9.+2.*math.sqrt(15))/21. ] ] )

def integrateFunctionOnUnitTriangle( f ):
  """This function returns the value of 

  \int_T f(x,y) d(x,y)

  where T is the right triangle with the cornes (0,0),(1,0),(0,1)."""
  return doIntegration( INTEGRATION_TRI_LOCATIONS, INTEGRATION_TRI_WEIGHTS, f )



  
INTEGRATION_LINE_LOCATIONS_3 = num.array( [
    0.1127016654,
    0.5,
    0.8872983346 ] )

INTEGRATION_LINE_WEIGHTS_3 = num.array( [
    5./18.,
    8./18.,
    5./18. ] )

def integrateFunctionOnUnitInterval3( f ):
  """This function returns the value of

  \int_0^1 f(x) dx.
  
  It uses a Gaussian quadrature of order 3."""
  return doIntegration( INTEGRATION_LINE_LOCATIONS_3, INTEGRATION_LINE_WEIGHTS_3, f )




INTEGRATION_LINE_LOCATIONS_4 = num.array( [
    0.0694318442,
    0.3300094782,
    0.6699905218,
    0.9305681558 ] )

INTEGRATION_LINE_WEIGHTS_4 = num.array( [
    0.1739274226,
    0.3260725774,
    0.3260725774,
    0.1739274226 ] )

def integrateFunctionOnUnitInterval4( f ):
  """This function returns the value of

  \int_0^1 f(x) dx.
  
  It uses a Gaussian quadrature of order 4."""
  return doIntegration( INTEGRATION_LINE_LOCATIONS_4, INTEGRATION_LINE_WEIGHTS_4, f )
