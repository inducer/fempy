import pylinear.matrices as num
import pylinear.linear_algebra as la
import math
import tools




class tErrorEstimator:
  def __init__(self, mesh, solution):
    self.Mesh = mesh
    self.Solution = solution
    self.ElementToErrorMap = {}

  def __call__(self, element):
    if element in self.ElementToErrorMap:
      return self.ElementToErrorMap[element]
    else:
      result = self._getEstimate(element)
      self.ElementToErrorMap[element] = result
      return result

  def estimateTotalError(self):
    return self.estimateErrorSum(self.Mesh.elements())

  def estimateErrorSum(self, elements):
    return sum([self(el) for el in elements])

  def _getEstimate(self, element):
    raise RuntimeError, "not implemented"




class tAnalyticSolutionL2ErrorEstimator(tErrorEstimator):
  def __init__(self, mesh, computed_solution, analytic_solution):
    tErrorEstimator.__init__(self, mesh, computed_solution)
    self.AnalyticSolution = analytic_solution

  def _getEstimate(self, element):
    def errorFunctionL2(point, solution_func_value):
      return (self.AnalyticSolution(element.transformToReal(point)) - solution_func_value) ** 2

    node_values = num.take(self.Solution, element.nodeNumbers())
    return element.getVolumeIntegralOver(errorFunctionL2, node_values)




class tAnalyticSolutionH1ErrorEstimator(tErrorEstimator):
  def __init__(self, mesh, computed_solution, analytic_solution, grad_analytic_solution):
    tErrorEstimator.__init__(self, mesh, computed_solution)
    self.AnalyticSolution = analytic_solution
    self.GradAnalyticSolution = grad_analytic_solution

    dims = mesh.dimensions()
    self.UnitVectors = []
    for i in range(dims):
      vec = num.zeros((dims,), num.Float)
      vec[i] = 1
      self.UnitVectors.append(vec)

  def _getEstimate(self, element):
    def errorFunctionH1(point, solution_func_value):
      real_point = element.transformToReal(point)
      xf_jac = element.getTransformJacobian(point)
      xf_jac_det = la.determinant(xf_jac)
      xf_jac_inv = la.inverse(xf_jac)
      twisted_grad = num.matrixmultiply(num.transpose(xf_jac_inv), approx_grad_solution(point))

      print "twisted:", twisted_grad
      print "exact:", self.GradAnalyticSolution(real_point)
      raw_input()

      return xf_jac_det * ((self.AnalyticSolution(real_point) - solution_func_value) ** 2 \
        + tools.norm2squared(twisted_grad-self.GradAnalyticSolution(real_point)))

    node_values = num.take(self.Solution, element.nodeNumbers())
    approx_grad_solution = element.getFormFunctionCombinationGradient(node_values)

    return element.getVolumeIntegralOver(errorFunctionH1, node_values)




