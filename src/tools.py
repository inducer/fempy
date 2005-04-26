import math, sys, operator, types

import pylinear.array as num
import pylinear.linear_algebra as la
import pylinear.operation as op
import pylinear.toybox as toybox




delta = toybox.delta




# Data structures ------------------------------------------------------------
class Reference(object):
    def __init__( self, value ):
        self.V = value
    def get( self ):
        return self.V
    def set( self, value ):
        self.V = value




class FunctionValueCache(object):
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

    


class LexicographicSequencer(object):
    def __init__(self, container, limits):
        self._Container = container
        self._Dimensions = [high-low for low, high in limits]
        self._Low = [low for low, high in limits]

    def __len__(self):
        return product(self._Dimensions)

    def translate_single_index(self, index):
        indices = []
        remaining_size = len(self)
        if not (0 <= index < remaining_size):
            raise IndexError, "invalid subscript to sequencer object"
        for i,sz in enumerate(self._Dimensions):
            remaining_size /= sz
            quotient, index = divmod(index, remaining_size)
            indices.append(quotient + self._Low[i])
        return tuple(indices)

    def get_all_indices(self):
        return [self.translate_single_index(i) for i in range(len(self))]

    def __getitem__(self, index):
        return self._Container[self.translate_single_index(index)]

  


class Grid(object):
    def __init__(self, origin, grid_vectors):
        self._Origin = origin
        self._GridVectors = grid_vectors

    def grid_vectors(self):
        return self._GridVectors

    def __getitem__(self, index):
        result = self._Origin.copy()
        for i, gv in zip(index, self._GridVectors):
            result += i * gv
        return result

    def find_closest_grid_point_index(self, point):
        tmat = num.array(self._GridVectors).T
        float_coords = tmat <<num.solve> (point - self._Origin)
        return tuple([int(round(c)) for c in float_coords])

    def interpolate_prid_point_index(self, point):
        tmat = num.array(self._GridVectors).T
        float_coords = tmat <<num.solve> (point - self._Origin)
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




class FiniteGrid(Grid):
    def __init__(self, origin, grid_vectors, limits):
        """Instantiates a finite grid. The limits are specified as a
        list of tuples of (low, high) values, one for each grid vector.
        For the index of a dimension, we assert, as is usual in Python:

        low <= index < high,

        such that there are (high-low) gridpoints and (high-low-1)
        grid intervals
        """
        assert len(grid_vectors) == len(limits)
        
        Grid.__init__(self, origin, grid_vectors)
        self._Limits = limits

    def limits(self):
        return self._Limits

    def __iter__(self):
        return iter(self.as_sequence().get_all_indices())

    def iterkeys(self):
        return self.__iter__()

    def grid_point_counts(self):
        """Returns the number of grid intervals in each direction.
        """
        return [high-low for low, high in self._Limits]

    def grid_point_count(self):
        """Returns the number of grid intervals in each direction.
        """
        return product(self.grid_point_counts())

    def is_within_bounds(self, key):
        for el, (low, high) in zip(key, self._Limits):
            if not (low <= el < high):
                return False
        return True

    def as_sequence(self):
        return LexicographicSequencer(self, self._Limits)

    def chop_upper_boundary(self, by = 1):
        return FiniteGrid(self._Origin, self._GridVectors,
                          [(low, high-by) for low, high in self._Limits])

    def chop_lower_boundary(self, by = 1):
        return FiniteGrid(self._Origin, self._GridVectors,
                          [(low+by, high) for low, high in self._Limits])

    def chop_both_boundaries(self, by = 1):
        return FiniteGrid(self._Origin, self._GridVectors,
                          [(low+by, high-by) for low, high in self._Limits])

    def enlarge_at_upper_boundary(self, by = 1):
        return FiniteGrid(self._Origin, self._GridVectors,
                          [(low, high+by) for low, high in self._Limits])

    def enlarge_at_lower_boundary(self, by = 1):
        return FiniteGrid(self._Origin, self._GridVectors,
                          [(low-by, high) for low, high in self._Limits])

    def enlarge_at_both_boundaries(self, by = 1):
        return FiniteGrid(self._Origin, self._GridVectors,
                          [(low-by, high+by) for low, high in self._Limits])

    def reduce_periodically(self, key):
        return tuple([
            el % (high-low) for el, (low, high) in zip(key, self._Limits)])

    def reduce_to_closest(self, key):
        return tuple([
            max(min(high-1, el), low) for el, (low, high) in zip(key, self._Limits)])
  



def make_subdivision_grid(origin, grid_vectors, limits):
    interval_counts = [high - low - 1 for low, high in limits]
    my_gvs = [gv / float(ivs) for gv, ivs in zip(grid_vectors, interval_counts)]
    return FiniteGrid(origin, my_gvs, limits)
    



def make_cell_centered_grid(origin, grid_vectors, limits):
    my_gvs = [gv / float(high - low) for gv, (low, high) in zip(grid_vectors, limits)]
    return FiniteGrid(origin + general_sum(my_gvs) * 0.5,
                       my_gvs, limits)
    



class DictionaryWithDefault(object):
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


    
class FakeList(object):
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




class DependentDictionary(object):
    def __init__(self, f, start = {}):
        self._Function = f
        self._Dictionary = start.copy()

    def copy(self):
        return DependentDictionary(self._Function, self._Dictionary)

    def __contains__(self, key):
        try:
            self[key]
            return True
        except KeyError:
            return False

    def __getitem__(self, key):
        try:
            return self._Dictionary[key]
        except KeyError:
            return self._Function(self._Dictionary, key)

    def __setitem__(self, key, value):
        self._Dictionary[key] = value
    
    def genuineKeys(self):
        return self._Dictionary.keys()

    def iteritems(self):
        return self._Dictionary.iteritems()

    def iterkeys(self):
        return self._Dictionary.iterkeys()

    def itervalues(self):
        return self._Dictionary.itervalues()




def add_tuples(t1, t2):
    return tuple([t1v + t2v for t1v, t2v in zip(t1, t2)])

def negate_tuple(t1):
    return tuple([-t1v for t1v in t1])




# Numerical algorithms -------------------------------------------------------
def abs_squared(x):
    return (x.conjugate() * x).real




def find_zero_by_newton(f, fprime, x_start, tolerance = 1e-12, maxit = 10):
    it = 0
    while it < maxit:
        it += 1
        f_value = f(x_start)
        if math.fabs(f_value) < tolerance:
            return x_start
        x_start -= f_value / fprime(x_start)
    raise RuntimeError, "Newton iteration failed, a zero was not found"




def find_vector_zero_by_newton(f, fprime, x_start, tolerance = 1e-12, maxit = 10):
    it = 0
    while it < maxit:
        it += 1
        f_value = f(x_start)
        if op.norm_2(f_value) < tolerance:
            return x_start
        x_start -= num.matrixmultiply(la.inverse(fprime(x_start)), f_value)
    raise RuntimeError, "Newton iteration failed, a zero was not found"




def distance_to_line(start_point, direction, point):
    # Ansatz: start_point + alpha * direction 
    # <start_point + alpha * direction - point, direction> = 0!
    alpha = - num.innerproduct(start_point - point, direction) / \
            op.norm_2_squared(direction)
    foot_point = start_point + alpha * direction
    return op.norm_2(point - foot_point), alpha




def angle_cosine_between_vectors(vec1, vec2):
    return vec1*vec2.H / (op.norm_2(vec1)*op.norm_2(vec2))




def interpolate_vector_list(vectors, inbetween_points):
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




def make_rotation_matrix(radians, n = 2, axis1 = 0, axis2 = 1, typecode = num.Float):
    mat = num.identity(n, typecode)
    mat[axis1,axis1] = math.cos(radians)
    mat[axis2,axis1] = math.sin(radians)
    mat[axis1,axis2] = -math.sin(radians)
    mat[axis2,axis2] = math.cos(radians)
    return mat




def get_parallelogram_volume(vectors):
    if vectors[0].shape[0] == 2:
        return vectors[0][0] * vectors[1][1] - vectors[1][0] * vectors[0][1]
    else:
        raise RuntimeError, "not implemented"




def unit_vector(i, dim, typecode = num.Float):
    uvec = num.zeros((dim,), typecode)
    uvec[i] = 1
    return uvec




def conjugate(value):
    try:
        return value.conjugate()
    except AttributeError:
        return value




class SparseVector(DictionaryWithDefault):
    def __init__(self):
        DictionaryWithDefault.__init__(self, lambda x: 0.)

    def add_to(self, other, factor = 1.):
        for key in self:
            other[key] += factor * self[key]

    def add_to_matrix_column(self, matrix, column, factor = 1.):
        for key in self:
            matrix[key, column] += factor * self[key]

    def add_to_matrix_row(self, matrix, row, factor = 1.):
        for key in self:
            matrix[row, key] += factor * self[key]

    def conjugate(self):
        result = SparseVector()
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
        result = SparseVector()
        for key in self:
            result[key] = other * self[key]
        return result

    def __rmul__(self, other):
        result = SparseVector()
        for key in self:
            result[key] = other * self[key]
        return result




class LinearSystemOfEquations:
    # UNTESTED.
    def __init__(self):
        self.Equations = []
        self.SymbolMap = {}
        pass

    def register_equation(self, coeffs_and_symbols, rhs):
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




def write_1d_gnuplot_graph(f, a, b, steps = 100, fname = ",,f.data", progress = False):
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

def write_1d_gnuplot_graphs(f, a, b, steps = 100, fnames = None, progress = False):
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



def write_gnuplot_graph(f, a, b, steps = 100, fname = ",,f.data", progress = False):
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




class DictionaryOfArithmeticTypes(dict):
    """Allows arithmetic operations on dictionaries
    which propagate to its elements.
    """

    def assert_same_keys(self, other):
        for key in self:
            assert key in other
        for key in other:
            assert key in self

    def unary_operator(self, operator):
        result = DictionaryOfArithmeticTypes()
        for key in self:
            result[key] = operator(self[key])
        return result

    def binary_operator(self, other, operator):
        try:
            self.assert_same_keys(other)
            result = DictionaryOfArithmeticTypes()
            for key in self:
                result[key] = operator(self[key], other[key])
            return result
        except TypeError:
            result = DictionaryOfArithmeticTypes()
            for key in self:
                result[key] = operator(self[key], other)
            return result

    def reverse_binary_operator(self, other, operator):
        try:
            self.assert_same_keys(other)
            result = DictionaryOfArithmeticTypes()
            for key in self:
                result[key] = operator(other[key], self[key])
            return result
        except TypeError:
            result = DictionaryOfArithmeticTypes()
            for key in self:
                result[key] = operator(other, self[key])
            return result

    def __neg__(self): return self.unary_operator(operator.neg)
    def __pos__(self): return self.unary_operator(operator.pos)
    def __abs__(self): return self.unary_operator(operator.abs)
    def __invert__(self): return self.unary_operator(operator.invert)

    def __add__(self, other): return self.binary_operator(other, operator.add)
    def __sub__(self, other): return self.binary_operator(other, operator.sub)
    def __mul__(self, other): return self.binary_operator(other, operator.mul)
    def __div__(self, other): return self.binary_operator(other, operator.div)
    def __mod__(self, other): return self.binary_operator(other, operator.mod)
    def __pow__(self, other): return self.binary_operator(other, operator.pow)
    def __lshift__(self, other): return self.binary_operator(other, operator.lshift)
    def __rshift__(self, other): return self.binary_operator(other, operator.rshift)
    def __and__(self, other): return self.binary_operator(other, operator.and_)
    def __or__(self, other): return self.binary_operator(other, operator.or_)
    def __xor__(self, other): return self.binary_operator(other, operator.xor)

    def __radd__(self, other): return self.reverse_binary_operator(other, operator.add)
    def __rsub__(self, other): return self.reverse_binary_operator(other, operator.sub)
    def __rmul__(self, other): return self.reverse_binary_operator(other, operator.mul)
    def __rdiv__(self, other): return self.reverse_binary_operator(other, operator.div)
    def __rmod__(self, other): return self.reverse_binary_operator(other, operator.mod)
    def __rpow__(self, other): return self.reverse_binary_operator(other, operator.pow)
    def __rlshift__(self, other): return self.reverse_binary_operator(other, operator.lshift)
    def __rrshift__(self, other): return self.reverse_binary_operator(other, operator.rshift)
    def __rand__(self, other): return self.reverse_binary_operator(other, operator.and_)
    def __ror__(self, other): return self.reverse_binary_operator(other, operator.or_)
    def __rxor__(self, other): return self.reverse_binary_operator(other, operator.xor)

    def __iadd__(self, other): 
        self.assert_same_keys(other)
        for key in self: 
            self[key] += other[key]
        return self

    def __isub__(self, other): 
        self.assert_same_keys(other)
        for key in self: 
            self[key] -= other[key]
        return self

    def __imul__(self, other): 
        self.assert_same_keys(other)
        for key in self: 
            self[key] *= other[key]
        return self

    def __idiv__(self, other): 
        self.assert_same_keys(other)
        for key in self: 
            self[key] /= other[key]
        return self

    def __imod__(self, other): 
        self.assert_same_keys(other)
        for key in self: 
            self[key] %= other[key]
        return self

    def __ipow__(self, other): 
        self.assert_same_keys(other)
        for key in self: 
            self[key] **= other[key]
        return self

    def __ilshift__(self, other): 
        self.assert_same_keys(other)
        for key in self: 
            self[key] <<= other[key]
        return self

    def __irshift__(self, other): 
        self.assert_same_keys(other)
        for key in self: 
            self[key] >>= other[key]
        return self

    def __iand__(self, other): 
        self.assert_same_keys(other)
        for key in self: 
            self[key] &= other[key]
        return self

    def __ior__(self, other): 
        self.assert_same_keys(other)
        for key in self: 
            self[key] |= other[key]
        return self

    def __ixor__(self, other): 
        self.assert_same_keys(other)
        for key in self: 
            self[key] ^= other[key]
        return self



# Generic utilities ----------------------------------------------------------
def flatten(list):
    result = []
    for i in list:
        result += i
    return result




def sum_over(function, arguments):
    result = 0
    for i in arguments:
        result += function(i)
    return result




def general_sum(sequence):
    return reduce(operator.add, sequence)




def linear_combination(coefficients, vectors):
    result = coefficients[0] * vectors[0]
    for c,v in zip(coefficients, vectors)[1:]:
        result += c*v
    return result




def average(sequence):
    return general_sum(sequence)/float(len(sequence))



def all_equal(sequence):
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




def cartesian_product(list1, list2):
    result = []
    for i in list1:
        for j in list2:
            result.append((i,j))




def cartesian_product_sum(list1, list2):
    """This routine returns a list of sums of each element of
    list1 with each element of list2. Also works with lists.
    """
    result = []
    for i in list1:
        for j in list2:
            result.append(i+j)
    return result




def reverse_dictionary(the_dict):
    result = {}
    for key, value in the_dict.iteritems():
        if value in result:
            raise RuntimeError, "non-reversible mapping"
        result[value] = key
    return result




def generate_positive_integer_tuples_below(n, length, least = 0):
    assert length >= 0
    if length == 0:
        yield []
    else:
        for i in range(least, n):
            for base in generate_positive_integer_tuples_below(n, length-1, least):
                yield [i] + base

def generate_all_positive_integer_tuples(length, least = 0):
    assert length >= 0
    current_max = least
    while True:
        for max_pos in range(length):
            for prebase in generate_positive_integer_tuples_below(current_max, max_pos, least):
                for postbase in generate_positive_integer_tuples_below(current_max+1, length-max_pos-1, least):
                    yield prebase + [current_max] + postbase
        current_max += 1

def _pos_and_neg_adaptor(tuple_iter):
    for tup in tuple_iter:
        nonzero_indices = [i for i in range(len(tup)) if tup[i] != 0]
        for do_neg_tup in generate_positive_integer_tuples_below(2, len(nonzero_indices)):
            this_result = list(tup)
            for index, do_neg in enumerate(do_neg_tup):
                if do_neg:
                    this_result[nonzero_indices[index]] *= -1
            yield tuple(this_result)

def generate_all_integer_tuples_below(n, length, least_abs = 0):
    return _pos_and_neg_adaptor(generate_positive_integer_tuples_below(
        n, length, least_abs))

def generate_all_integer_tuples(length, least_abs = 0):
    return _pos_and_neg_adaptor(generate_all_positive_integer_tuples(
        length, least_abs))
            



# Obscure stuff --------------------------------------------------------------
def write_matrix_as_csv(filename, matrix):
    mat_file = file(filename, "w")
    h,w = matrix.shape
    for row in range(0, h):
        for column in range(0, w):
            mat_file.write("%f," % matrix[ row, column ])
    mat_file.write("\n")




def enumerate_basic_directions(dimensions):
    coordinate_list = [[0], [1], [-1]]
    return reduce(cartesian_product_sum, [coordinate_list] * dimensions)[1:]



