import weakref
import pylinear.matrices as num
import pylinear.matrix_tools as mtools
import pylinear.linear_algebra as la
import stopwatch
import tools




def getValueOnElement(element, node_v, unit_point):
    ffs = element.formFunctionKit().formFunctions()
    nonus = element.nodeNumbers()

    value = ffs[0](unit_point) * node_v[nonus[0]]

    for ff, nonu in zip(ffs, nonus)[1:]:
        value += ff(unit_point) * node_v[nonu]
    return value




def getGradientOnElement(element, node_v, unit_point):
    ffs = element.formFunctionKit().differentiatedFormFunctions()
    nonus = element.nodeNumbers()

    value = num.array([deriv(unit_point) for deriv in ffs[0]]) * node_v[nonus[0]]

    for grad, nonu in zip(ffs, nonus)[1:]:
        value += num.array([deriv(unit_point) for deriv in grad]) * node_v[nonu]
    return value




def getRealGradientOnElement(element, node_v, unit_point):
    return num.matrixmultiply(
        num.transpose(la.inverse(
        element.getTransformJacobian(unit_point))),
        getGradientOnElement(element, node_v, unit_point))




def getGlobalSolution(mesh, node_v, point):
    el = mesh.findElement(point)
    return getValueOnElement(el, node_v, el.transformToUnit(point))




def getGlobalSolutionGradient(mesh, node_v, point):
    el = mesh.findElement(point)
    return getRealGradientOnElement(el, node_v, 
                                    el.transformToUnit(point))




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
