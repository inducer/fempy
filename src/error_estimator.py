import pylinear.matrices as num
import math




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
      return (self.AnalyticSolution(point) - solution_func_value) ** 2

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
    def approxFirstDerivative(f, point, direction):
      # centered difference
      return (f(point + scale * direction) - f(point - scale * direction)) / (2 * scale)

    # This is incredibly slow. But who cares.
    approx_solution = element.getSolutionFunction(self.Solution)
    scale = math.sqrt(element.area()) * 0.01

    def errorFunctionH1(point, solution_func_value):
      approx_solution_diff = [approxFirstDerivative(approx_solution, point, vec) 
        for vec in self.UnitVectors]
      return (self.AnalyticSolution(point) - solution_func_value) ** 2 \
        + sum([(ana_diff(point) - approx_diff)**2 for ana_diff, approx_diff in 
            zip(self.GradAnalyticSolution, approx_solution_diff)])

    node_values = num.take(self.Solution, element.nodeNumbers())
    return element.getVolumeIntegralOver(errorFunctionH1, node_values)




