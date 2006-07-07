import pylinear.array as num
import pylinear.operation as op
import math
import pytools

# Weights from H.R. Schwarz, "Methode der finiten Elemente", Teubner 1980,
# p. 121f.

_INTEGRATION_TRI_WEIGHTS = num.array([
    0.1125, # barycenter
    
    (155.+math.sqrt(15))/2400.,
    (155.+math.sqrt(15))/2400.,
    (155.+math.sqrt(15))/2400.,

    (155.-math.sqrt(15))/2400.,
    (155.-math.sqrt(15))/2400.,
    (155.-math.sqrt(15))/2400. ])

_INTEGRATION_TRI_LOCATIONS = [
    num.array([ 1/3.,1/3., ]),

    num.array([ (6.+math.sqrt(15))/21., (6.+math.sqrt(15))/21. ]),
    num.array([ (9.-2.*math.sqrt(15))/21., (6.+math.sqrt(15))/21. ]),
    num.array([ (6.+math.sqrt(15))/21., (9.-2.*math.sqrt(15))/21. ]),

    num.array([ (6.-math.sqrt(15))/21., (6.-math.sqrt(15))/21. ]),
    num.array([ (9.+2.*math.sqrt(15))/21., (6.-math.sqrt(15))/21. ]),
    num.array([ (6.-math.sqrt(15))/21., (9.+2.*math.sqrt(15))/21. ]) ]

_INTEGRATION_TRI = zip(_INTEGRATION_TRI_WEIGHTS,
                       _INTEGRATION_TRI_LOCATIONS)

def integrate_on_unit_triangle(f):
    """This function returns the value of 
    
    \int_T f(x,y) d(x,y)
    
    where T is the right triangle with the cornes (0,0),(1,0),(0,1)."""
    result = 0
    for weight, location in _INTEGRATION_TRI:
        result += weight * f(location)
    return result

def integrate_on_unit_triangle_with_order_1(f):
    return 0.5 * f(num.array([1/3., 1/3.]))




  
_INTEGRATION_LINE_LOCATIONS_3 = num.array([
    0.1127016654,
    0.5,
    0.8872983346 ])

_INTEGRATION_LINE_WEIGHTS_3 = num.array([
    5./18.,
    8./18.,
    5./18. ])

_INTEGRATION_LINE_3 = zip(_INTEGRATION_LINE_WEIGHTS_3,
                          _INTEGRATION_LINE_LOCATIONS_3)

def integrate_on_unit_interval_3(f):
    """This function returns the value of
    
    \int_0^1 f(x) dx.
    
    It uses a Gaussian quadrature of order 3."""
    result = 0
    for weight, location in _INTEGRATION_LINE_3:
        result += weight * f(location)
    return result




_INTEGRATION_LINE_LOCATIONS_4 = num.array([
    0.0694318442,
    0.3300094782,
    0.6699905218,
    0.9305681558 ])

_INTEGRATION_LINE_WEIGHTS_4 = num.array([
    0.1739274226,
    0.3260725774,
    0.3260725774,
    0.1739274226 ])

_INTEGRATION_LINE_4 = zip(_INTEGRATION_LINE_WEIGHTS_4,
                          _INTEGRATION_LINE_LOCATIONS_4)

def integrate_on_unit_interval_4(f):
    """This function returns the value of
    
    \int_0^1 f(x) dx.
    
    It uses a Gaussian quadrature of order 4."""
    result = 0
    for weight, location in _INTEGRATION_LINE_4:
        result += weight * f(location)
    return result




def integrate_along_line(point1, point2, f):
    direction = point2 - point1
    return op.norm_2(direction)\
           * integrate_on_unit_interval_4(lambda x: f(point1 + x * direction))




def integrate_on_two_dimensional_grid(grid, f):
    dims = grid.chop_upper_boundary().grid_point_counts()
    assert len(dims) == 2
    dim1 = dims[0]
    dim2 = dims[1]
    result = .25 * (f((0,0)) + f((dim1,0)) + f((0,dim2)) + f((dim1,dim2)))
    for i in range(1,dim1):
        for j in range(1,dim2):
            result += f((i,j))
    for i in range(1,dim1):
        result += .5 * (f((i,0))+f((i,dim2)))
    for j in range(1,dim2):
        result += .5 * (f((0,j))+f((dim1,j)))

    gv = grid.grid_vectors()
    return result * pytools.get_parallelogram_volume(gv)
