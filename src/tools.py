import math, sys, operator, types

import pylinear.matrices as num
import pylinear.linear_algebra as la
import pylinear.matrix_tools as mtools




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

    def iterkeys(self):
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

    def chopUpperBoundary(self):
        return tFiniteGrid(self._Origin, self._GridVectors,
                           [(low, high-1) for low, high in self._Limits])

    def chopLowerBoundary(self):
        return tFiniteGrid(self._Origin, self._GridVectors,
                           [(low+1, high) for low, high in self._Limits])

    def chopBothBoundaries(self):
        return tFiniteGrid(self._Origin, self._GridVectors,
                           [(low+1, high-1) for low, high in self._Limits])

    def reducePeriodically(self, key):
        return tuple([
            el % (high-low) for el, (low, high) in zip(key, self._Limits)])
  



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
def absSquared(x):
    return (x.conjugate() * x).real




def findZeroByNewton(f, fprime, x_start, tolerance = 1e-12, maxit = 10):
    it = 0
    while it < maxit:
        it += 1
        f_value = f(x_start)
        if math.fabs(f_value) < tolerance:
            return x_start
        x_start -= f_value / fprime(x_start)
    raise RuntimeError, "Newton iteration failed, a zero was not found"




def findVectorZeroByNewton(f, fprime, x_start, tolerance = 1e-12, maxit = 10):
    it = 0
    while it < maxit:
        it += 1
        f_value = f(x_start)
        if norm2(f_value) < tolerance:
            return x_start
        x_start -= num.matrixmultiply(la.inverse(fprime(x_start)), f_value)
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
    uvec = num.zeros((dim,), typecode)
    uvec[i] = 1
    return uvec




def conjugate(value):
    try:
        return value.conjugate()
    except AttributeError:
        return value




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




class tLinearSystemOfEquations:
    # UNTESTED.
    def __init__(self):
        self.Equations = []
        self.SymbolMap = {}
        pass

    def registerEquation(self, coeffs_and_symbols, rhs):
        self.Equations.append(([
            (coeff, self.SymbolMap.setdefault(symbol, len(self.SymbolMap)))
            for coeff, symbol  in coeffs_and_symbols], rhs))

    def solve(self, typecode = num.Float):
        m = num.zeros((len(self.Equations), len(self.SymbolMap)), typecode)
        rhs = num.zeros((len(self.SymbolMap),), typecode)
        for i, (eq, rhs) in enumerate(self.Equations):
            for coeff, j in eq:
                m[i,j] = coeff
            rhs[i] = rhs
        sol = la.solve_linear_equations(m, rhs)
        result = {}
        for sym, index in self.SymbolMap.iteritems():
            result[sym] = sol[index]
        return result




def write1DGnuplotGraph(f, a, b, steps = 100, fname = ",,f.data", progress = False):
    h = float(b - a)/steps
    gnuplot_file = file(fname, "w")

    def do_plot(func):
        for n in range(steps):
            if progress:
                sys.stdout.write(".")
                sys.stdout.flush()
            x = a + h * n
            gnuplot_file.write("%f\t%f\n" % (x, func(x)))

    do_plot(f)
    if progress:
        sys.stdout.write("\n")

def write1DGnuplotGraphs(f, a, b, steps = 100, fnames = None, progress = False):
    h = float(b - a)/steps
    if not fnames:
        result_count = len(f(a))
        fnames = [",,f%d.data" % i for i in range(result_count)]

    gnuplot_files = [file(fname, "w") for fname in fnames]

    for n in range(steps):
        if progress:
            sys.stdout.write(".")
            sys.stdout.flush()
        x = a + h * n
        for gpfile, y in zip(gnuplot_files, f(x)):
            gpfile.write("%f\t%f\n" % (x, y))
    if progress:
        sys.stdout.write("\n")



def writeGnuplotGraph(f, a, b, steps = 100, fname = ",,f.data", progress = False):
    h = float(b - a)/steps
    gnuplot_file = file(fname, "w")

    def do_plot(func):
        for n in range(steps):
            if progress:
                sys.stdout.write(".")
                sys.stdout.flush()
            x = a + h * n
            gnuplot_file.write("%f\t%f\n" % (x, func(x)))

    if isinstance(f, types.ListType):
        for f_index, real_f in enumerate(f):
            if progress:
                sys.stdout.write("function %d: " % f_index)
            do_plot(real_f)
            gnuplot_file.write("\n")
            if progress:
                sys.stdout.write("\n")
    else:
        do_plot(f)
        if progress:
            sys.stdout.write("\n")




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




def generalSum(sequence):
    return reduce(operator.add, sequence)




def linearCombination(coefficients, vectors):
    result = coefficients[0] * vectors[0]
    for c,v in zip(coefficients, vectors)[1:]:
        result += c*v
    return result




def average(sequence):
    return generalSum(sequence)/float(len(sequence))



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




def argmax(list, f = lambda x: x):
    current_max_index = -1
    current_max = f(list[0])

    for idx, item in enumerate(list[1:]):
        value = f(item)
        if value > current_max:
            current_max_index = idx
            current_max = value
    return current_max_index+1




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




def generateIntegerTuplesBelow(n, length, least = 0):
    assert length >= 0
    if length == 0:
        yield []
    else:
        for i in range(least, n):
            for base in generateIntegerTuplesBelow(n, length-1, least):
                yield [i] + base

def generateAllIntegerTuples(length, least = 0):
    assert length >= 0
    current_max = least
    while True:
        for max_pos in range(length):
            for prebase in generateIntegerTuplesBelow(current_max, max_pos, least):
                for postbase in generateIntegerTuplesBelow(current_max+1, length-max_pos-1, least):
                    yield prebase + [current_max] + postbase
        current_max += 1
            



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



