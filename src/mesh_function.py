import weakref
import element
import pylinear.matrices as num
import pylinear.matrix_tools as mtools
import pylinear.linear_algebra as la
import stopwatch
import tools




class tMeshFunction(object):
    def __init__(self, mesh, number_assignment, vector):
        self.Mesh = mesh
        self.NumberAssignment = number_assignment
        self.Vector = vector

    # accessors ------------------------------------------------------------------
    def mesh(self):
        return self.Mesh

    def numberAssignment(self):
        return self.NumberAssignment

    def vector(self):
        return self.Vector

    # value getters --------------------------------------------------------------
    def __getitem__(self, node):
        return element.getNodeValue(node, self.NumberAssignment, self.Vector)

    def getValueOnElement(self, el, unit_point):
        ffs = el.formFunctionKit().formFunctions()
        nodes = el.nodes()

        value = ffs[0](unit_point) * element.getNodeValue(nodes[0], 
                                                          self.NumberAssignment,
                                                          self.Vector)
        
        for ff, node in zip(ffs, nodes)[1:]:
            value += ff(unit_point) * element.getNodeValue(node,
                                                           self.NumberAssignment,
                                                           self.Vector)
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
