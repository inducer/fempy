import weakref
import element
import pylinear.matrices as num
import pylinear.matrix_tools as mtools
import pylinear.linear_algebra as la
import stopwatch
import tools




class tMeshFunction(object):
    def __init__(self, mesh, number_assignment, vector, node_constraints):
        self.Mesh = mesh
        self.NumberAssignment = number_assignment
        self.Vector = vector
        self.NodeConstraints = node_constraints

    def copy(self, mesh = None, number_assignment = None, vector = None,
             node_constraints = None):
        return tMeshFunction(
            mesh or self.Mesh,
            number_assignment or self.NumberAssignment,
            vector or self.Vector,
            node_constraints or self.NodeConstraints)

    # accessors ------------------------------------------------------------------
    def mesh(self):
        return self.Mesh

    def numberAssignment(self):
        return self.NumberAssignment

    def vector(self):
        return self.Vector

    # operations -----------------------------------------------------------------
    def conjugate(self):
        return tMeshFunction(self.Mesh, self.NumberAssignment, 
                             num.conjugate(self.Vector))

    def __add__(self, other):
        assert self.Mesh is other.Mesh
        assert self.NumberAssignment is other.NumberAssignment
        return tMeshFunction(self.Mesh, self.NumberAssignment, 
                             self.Vector + other.Vector)

    def __sub__(self, other):
        assert self.Mesh is other.Mesh
        assert self.NumberAssignment is other.NumberAssignment
        return tMeshFunction(self.Mesh, self.NumberAssignment, 
                             self.Vector - other.Vector)

    def __mul__(self, factor):
        return self.copy(vector = factor * self.Vector)

    def __rmul__(self, factor):
        return self.copy(vector = factor * self.Vector)

    def realPart(self):
        return self.copy(vector = self.Vector.real)

    def imaginaryPart(self):
        return self.copy(vector = self.Vector.imaginary)

    real = property(realPart)
    imaginary = property(imaginaryPart)

    # value getters --------------------------------------------------------------
    def __getitem__(self, node):
        try:
            constraint = self.NodeConstraints[node]
            return constraint[0] \
                   + sum([coeff * getNodeValue(other_node, number_assignment, vector)
                          for coeff, other_node in constraint[1]])
        except KeyError:
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

        value = num.array([deriv(unit_point) for deriv in ffs[0]]) * self[nodes[0]]
                               
        for grad, node in zip(ffs, nodes)[1:]:
            value += num.array([deriv(unit_point) for deriv in grad]) * self[node]
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




class tScalarProductCalculator:
    def __init__(self, mass_matrix):
        self.MassMatrix = mass_matrix
        self.CastMassMatrix = {}

    def __getinitargs__(self):
        return self.MassMatrix,

    def __getstate__(self):
        return None

    def __setstate__(self, state):
        pass
    
    def __call__(self, v1, v2):
        tc = v1.typecode()
        try:
            cast_spi = self.CastMassMatrix[tc]
        except:
            cast_spi = num.asarray(self.MassMatrix, tc, num.SparseExecuteMatrix)
            self.CastMassMatrix[tc] = cast_spi

        return mtools.sp(v1, num.matrixmultiply(cast_spi, v2))
