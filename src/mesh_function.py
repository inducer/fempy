import pylinear.array as num
import pylinear.linear_algebra as la




def makeCompleteVector(complete_number_assignment, incomplete_vector,
                       constraints):
    assert min(complete_number_assignment.values()) == 0
    assert max(complete_number_assignment.values()) == len(complete_number_assignment) - 1

    incomplete_dof_count = incomplete_vector.shape[0]
    complete_dof_count = len(complete_number_assignment)

    result = num.zeros((complete_dof_count,), incomplete_vector.typecode())
    result[:incomplete_dof_count] = incomplete_vector
    
    def get_value(node):
        try:
            constraint = constraints[node]
            return constraint[0] \
                   + sum([coeff * get_value(other_node)
                          for coeff, other_node in constraint[1]])
        except KeyError:
            return incomplete_vector[complete_number_assignment[node]]

    for node in constraints:
        result[complete_number_assignment[node]] = get_value(node)
    return result




class tMeshFunction(object):
    def __init__(self, mesh, number_assignment, vector):
        self.Mesh = mesh
        self.NumberAssignment = number_assignment
        self.Vector = vector

    def copy(self, mesh = None, number_assignment = None, vector = None):
        return tMeshFunction(
            mesh or self.Mesh,
            number_assignment or self.NumberAssignment,
            vector or self.Vector)

    # accessors ------------------------------------------------------------------
    def mesh(self):
        return self.Mesh

    def numberAssignment(self):
        return self.NumberAssignment

    def vector(self):
        return self.Vector

    # operations -----------------------------------------------------------------
    def conjugate(self):
        return self.copy(vector = num.conjugate(self.Vector))

    def __add__(self, other):
        assert self.Mesh is other.Mesh
        assert self.NumberAssignment is other.NumberAssignment

        return self.copy(vector = self.Vector + other.Vector)

    def __iadd__(self, other):
        assert self.Mesh is other.Mesh
        assert self.NumberAssignment is other.NumberAssignment

        self.Vector += other.Vector
        return self

    def __sub__(self, other):
        assert self.Mesh is other.Mesh
        assert self.NumberAssignment is other.NumberAssignment

        return self.copy(vector = self.Vector - other.Vector)

    def __isub__(self, other):
        assert self.Mesh is other.Mesh
        assert self.NumberAssignment is other.NumberAssignment

        self.Vector -= other.Vector
        return self

    def __mul__(self, factor):
        return self.copy(vector = factor * self.Vector)

    def __imul__(self, factor):
        self.Vector *= factor
        return self

    def __rmul__(self, factor):
        return self.copy(vector = factor * self.Vector)

    def __div__(self, divisor):
        return self.copy(vector = self.Vector/ divisor)

    def __idiv__(self, divisor):
        self.Vector /= divisor
        return self

    def realPart(self):
        return self.copy(vector = self.Vector.real)

    def imaginaryPart(self):
        return self.copy(vector = self.Vector.imaginary)

    real = property(realPart)
    imaginary = property(imaginaryPart)

    # value getters --------------------------------------------------------------
    def __getitem__(self, node):
        return self.Vector[self.NumberAssignment[node]]

    def getValueOnElement(self, el, unit_point):
        ffs = el.formFunctionKit().formFunctions()
        nodes = el.nodes()

        value = ffs[0](unit_point) * self[nodes[0]]
        
        for ff, node in zip(ffs, nodes)[1:]:
            value += ff(unit_point) * self[node]
        return value

    def getGradientOnElement(self, el, unit_point):
        ffs = el.formFunctionKit().differentiatedFormFunctions()
        nodes = el.nodes()

        value =  num.array([self[nodes[0]] * deriv(unit_point) for deriv in ffs[0]])
                               
        for grad, node in zip(ffs, nodes)[1:]:
            value += num.array([self[node] * deriv(unit_point) for deriv in grad])
        return value

    def getRealGradientOnElement(self, el, unit_point):
        return num.matrixmultiply(
            num.transpose(la.inverse(
            el.getTransformJacobian(unit_point))),
            self.getGradientOnElement(el, unit_point))

    def __call__(self, point):
        el = self.Mesh.findElement(point)
        return self.getValueOnElement(el, el.transformToUnit(point))

    def getGradient(self, point):
        el = self.Mesh.findElement(point)
        return self.getRealGradientOnElement(el, el.transformToUnit(point))




def discretizeFunction(mesh, f, typecode = num.Float, number_assignment = None):
    if number_assignment is None:
        number_assignment = {}
        for i, node in enumerate(mesh.dofManager()):
            number_assignment[node] = i

    vector = num.zeros((len(mesh.dofManager()),), typecode)
    for node in mesh.dofManager():
        vector[number_assignment[node]] = f(node.Coordinates)

    return tMeshFunction(mesh, number_assignment, vector)




def discretizeFunctionByIntegrals(mesh, f, mass_matrix, number_assignment, 
                       typecode = num.Float):
    # use the formula:
    # \sum_j\alpha_j \int\psi_i(x)\psi_j(x) dx =!= \int f(x) \psi_i(x)
    # <=>
    # M \alpha = \beta
    # where \alpha and \beta are vectors.

    mass_matrix = num.asarray(mass_matrix, typecode)
    beta = num.zeros((len(mesh.dofManager()),), typecode)
    for el in mesh.elements():
        vol_ints = el.getVolumeIntegralsOverFormFunction(lambda x, phi_x: f(x)*phi_x,
                                                         typecode)
        for node, vol_int in zip(el.nodes(), vol_ints):
            beta[number_assignment[node]] += vol_int

    alpha = mass_matrix <<num.solve>> beta
    
    return tMeshFunction(mesh, number_assignment, alpha)




class tScalarProductCalculator:
    def __init__(self, number_assignment, mass_matrix):
        self.MassMatrix = mass_matrix
        self.NumberAssignment = number_assignment
        self.CastMassMatrix = {}

    def __getinitargs__(self):
        return self.MassMatrix,

    def __getstate__(self):
        return None

    def __setstate__(self, state):
        pass
    
    def __call__(self, mf1, mf2):
        v1 = mf1.vector()
        v2 = mf2.vector()

        assert mf1.numberAssignment() is self.NumberAssignment
        assert mf2.numberAssignment() is self.NumberAssignment

        tc = v1.typecode()
        try:
            cast_spi = self.CastMassMatrix[tc]
        except:
            cast_spi = num.asarray(self.MassMatrix, tc, num.SparseExecuteMatrix)
            self.CastMassMatrix[tc] = cast_spi

        return v1 * (cast_spi * v2.H)
