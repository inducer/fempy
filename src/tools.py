import pylinear.matrices as num
import pylinear.linear_algebra as la
import pylinear.matrix_tools as mtools
import math




norm2 = mtools.norm2
norm2squared = mtools.norm2squared
delta = mtools.delta




# Data structures ------------------------------------------------------------
class tReference(object):
    def __init__( self, value ):
        self.V = value
    def get( self ):
        return self.V
    def set( self, value ):
        self.V = value




class tFunctionValueCache(object):
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

    


class tLexicographicSequencer(object):
    def __init__(self, container, limits):
        self._Container = container
        self._Dimensions = [high-low for low, high in limits]
        self._Low = [low for low, high in limits]

    def __len__(self):
        return product(self._Dimensions)

    def translateSingleIndex(self, index):
        indices = []
        remaining_size = len(self)
        if not (0 <= index < remaining_size):
            raise IndexError, "invalid subscript to sequencer object"
        for i,sz in enumerate(self._Dimensions):
            remaining_size /= sz
            quotient, index = divmod(index, remaining_size)
            indices.append(quotient + self._Low[i])
        return tuple(indices)

    def getAllIndices(self):
        return [self.translateSingleIndex(i) for i in range(len(self))]

    def __getitem__(self, index):
        return self._Container[self.translateSingleIndex(index)]

  


class tGrid(object):
    def __init__(self, origin, grid_vectors):
        self._Origin = origin
        self._GridVectors = grid_vectors

    def gridVectors(self):
        return self._GridVectors

    def __getitem__(self, index):
        result = self._Origin.copy()
        for i, gv in zip(index, self._GridVectors):
            result += i * gv
        return result

    def findClosestGridPointIndex(self, point):
        tmat = num.transpose(num.array(self._GridVectors))
        float_coords = la.solve_linear_equations(tmat, point - self._Origin)
        return tuple([int(round(c)) for c in float_coords])

    def interpolateGridPointIndex(self, point):
        tmat = num.transpose(num.array(self._GridVectors))
        float_coords = la.solve_linear_equations(tmat, point - self._Origin)
        rounded_down_int_coords = [int(math.floor(c)) for c in float_coords]
        neighbors = [rounded_down_int_coords]
        for d in range(len(self._GridVectors)):
            new_neighbors = []
            for item in neighbors:
                new_neighbors.append(item) 
                new_neighbor = item[:]
                new_neighbor[d] += 1
                new_neighbors.append(new_neighbor)
            neighbors = new_neighbors
        weights = []
        for neighbor in neighbors:
            weight = product([1-abs(a-b) for a,b in zip(float_coords, neighbor)])
            if abs(weight) >= 1e-5:
                weights.append((weight, tuple(neighbor)))
        return weights




class tFiniteGrid(tGrid):
    def __init__(self, origin, grid_vectors, limits):
        """Instantiates a finite grid. The limits are specified as a
        list of tuples of (low, high) values, one for each grid vector.
        For the index of a dimension, we assert, as is usual in Python:

        low <= index < high,

        such that there are (high-low) gridpoints and (high-low-1)
        grid intervals
        """
        assert len(grid_vectors) == len(limits)
        
        tGrid.__init__(self, origin, grid_vectors)
        self._Limits = limits

    def limits(self):
        return self._Limits

    def __iter__(self):
        return iter(self.asSequence().getAllIndices())

    def iterkeys():
        return self.__iter__()

    def gridPointCounts(self):
        """Returns the number of grid intervals in each direction.
        """
        return [high-low for low, high in self._Limits]

    def gridIntervalCounts(self):
        """Returns the number of grid intervals in each direction.
        """
        return [high-low-1 for low, high in self._Limits]
    
    def asSequence(self):
        return tLexicographicSequencer(self, self._Limits)

    def gridBlockIndices(self):
        seq = tLexicographicSequencer(self, [(low, high-1) for low, high in self._Limits])
        for idx in range(len(seq)):
            multidim_idx = seq.translateSingleIndex(idx)
            yield multidim_idx 
  



def makeSubdivisionGrid(origin, grid_vectors, limits):
    interval_counts = [high - low - 1 for low, high in limits]
    my_gvs = [gv / float(ivs) for gv, ivs in zip(grid_vectors, interval_counts)]
    return tFiniteGrid(origin, my_gvs, limits)
    



class tDictionaryWithDefault(object):
    def __init__(self, default_value_generator, start = {}):
        self._Dictionary = dict(start)
        self._DefaultGenerator = default_value_generator

    def __getitem__(self, index):
        try:
            return self._Dictionary[index]
        except KeyError:
            value = self._DefaultGenerator(index)
            self._Dictionary[index] = value
            return value

    def __setitem__(self, index, value):
        self._Dictionary[index] = value

    def __contains__(self, item):
        return True

    def iterkeys(self):
        return self._Dictionary.iterkeys()

    def __iter__(self):
        return self._Dictionary.__iter__()

    def iteritems(self):
        return self._Dictionary.iteritems()


    
class tFakeList(object):
    def __init__(self, f, length):
        self._Length = length
        self._Function = f

    def __len__(self):
        return self._Length

    def __getitem__(self, index):
        try:
            return [self._Function(i)
                    for i in range(*index.indices(self._Length))]
        except AttributeError:
            return self._Function(index)




class tDependentDictionary(object):
    def __init__(self, f, start = {}):
        self._Function = f
        self._Dictionary = start.copy()

    def copy(self):
        return tDependentDictionary(self._Function, self._Dictionary)

    def __getitem__(self, key):
        try:
            return self._Dictionary[key]
        except KeyError:
            return self._Function(self._Dictionary, key)

    def __setitem__(self, key, value):
        self._Dictionary[key] = value
    
    def genuineKeys(self):
        return self._Dictionary.keys()




# Numerical algorithms -------------------------------------------------------
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




def distanceToLine(start_point, direction, point):
    # Ansatz: start_point + alpha * direction 
    # <start_point + alpha * direction - point, direction> = 0!
    alpha = - num.innerproduct(start_point - point, direction) / \
            norm2squared(direction)
    foot_point = start_point + alpha * direction
    return norm2(point - foot_point), alpha




def angleCosineBetweenVectors(vec1, vec2):
    return mtools.sp(vec1, vec2) / (norm2(vec1)*norm2(vec2))




def interpolateVectorList(vectors, inbetween_points):
    if len(vectors) == 0:
        return []

    result = [vectors[0]]
    last_vector = vectors[0]
    for vector in vectors[1:]:
        for i in range(inbetween_points):
            result.append(last_vector + (vector-last_vector) \
                          * float(i+1) \
                          / float(inbetween_points+1))
        result.append(vector)
        last_vector = vector
    return result




def makeRotationMatrix(radians, n = 2, axis1 = 0, axis2 = 1, typecode = num.Float):
    mat = num.identity(n, typecode)
    mat[axis1,axis1] = math.cos(radians)
    mat[axis2,axis1] = math.sin(radians)
    mat[axis1,axis2] = -math.sin(radians)
    mat[axis2,axis2] = math.cos(radians)
    return mat




def getParallelogramVolume(vectors):
    if vectors[0].shape[0] == 2:
        return vectors[0][0] * vectors[1][1] - vectors[1][0] * vectors[0][1]
    else:
        raise RuntimeError, "not implemented"




def unitVector(i, dim, typecode = num.Float):
    uvec = num.zeros((dim,), num.Float)
    uvec[i] = 1
    return uvec




def conjugate(value):
    try:
        return value.conjugate()
    except AttributeError:
        return value




def frobeniusNorm(a):
    result = 0
    for i,j in a.indices():
        result += abs(a[i,j])**2
    return result

def matrixExp(a, eps = 1e-15):
    h,w = a.shape
    assert h == w
    a_frob = frobeniusNorm(a)
    
    last_result = num.identity(h, a.typecode())
    result = last_result.copy()

    current_power_of_a = a

    factorial = 1
    n = 1

    while True:
        result += current_power_of_a * (1./factorial)

        if frobeniusNorm(result - last_result)/a_frob < eps:
            return result

        n += 1
        last_result = result.copy()
        factorial *= n
        current_power_of_a = num.matrixmultiply(current_power_of_a, a)
    
        
    

class tSparseVector(tDictionaryWithDefault):
    def __init__(self):
        tDictionaryWithDefault.__init__(self, lambda x: 0.)

    def addTo(self, other, factor = 1.):
        for key in self:
            other[key] += factor * self[key]

    def addToMatrixColumn(self, matrix, column, factor = 1.):
        for key in self:
            matrix[key, column] += factor * self[key]

    def addToMatrixRow(self, matrix, row, factor = 1.):
        for key in self:
            matrix[row, key] += factor * self[key]

    def conjugate(self):
        result = tSparseVector()
        for key in self:
            result[key] = (self[key]+0j).conjugate()
        return result

    def __radd__(self, other):
        result = other.copy()
        for key in self:
            result[key] += self[key]
        return result

    def __add__(self, other):
        result = other.copy()
        for key in self:
            result[key] += self[key]
        return result

    def __rsub__(self, other):
        result = other.copy()
        for key in self:
            result[key] -= self[key]
        return result

    def __sub__(self, other):
        result = other.copy()
        for key in self:
            result[key] = self[key] - result[key]
        return result

    def __mul__(self, other):
        result = tSparseVector()
        for key in self:
            result[key] = other * self[key]
        return result

    def __rmul__(self, other):
        result = tSparseVector()
        for key in self:
            result[key] = other * self[key]
        return result




# Generic utilities ----------------------------------------------------------
def flatten(list):
    result = []
    for i in list:
        result += i
    return result




def sumOver(function, arguments):
    result = 0
    for i in arguments:
        result += function(i)
    return result




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




def product(list):
    return reduce(lambda x,y: x*y, list, 1)




def argmin(list, f = lambda x: x):
    current_min_index = -1
    current_min = f(list[0])

    for idx, item in enumerate(list[1:]):
        value = f(item)
        if value < current_min:
            current_min_index = idx
            current_min = value
    return current_min_index+1




def cartesianProduct(list1, list2):
    result = []
    for i in list1:
        for j in list2:
            result.append((i,j))




def cartesianProductSum(list1, list2):
    """This routine returns a list of sums of each element of
    list1 with each element of list2. Also works with lists.
    """
    result = []
    for i in list1:
        for j in list2:
            result.append(i+j)
    return result




def reverseDictionary(the_dict):
    result = {}
    for key, value in the_dict.iteritems():
        if value in result:
            raise RuntimeError, "non-reversible mapping"
        result[value] = key
    return result




# Obscure stuff --------------------------------------------------------------
def writeMatrixAsCSV(filename, matrix):
    mat_file = file(filename, "w")
    h,w = matrix.shape
    for row in range(0, h):
        for column in range(0, w):
            mat_file.write("%f," % matrix[ row, column ])
    mat_file.write("\n")




def enumerateBasicDirections(dimensions):
    coordinate_list = [[0], [1], [-1]]
    return reduce(cartesianProductSum, [coordinate_list] * dimensions)[1:]




