import pylinear.matrices as num




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

  def estimateTotalError():
    return self.estimateErrorSum(self.Mesh.elements())

  def estimateErrorSum(elements):
    return sum([self.estimateError(el) for el in elements])

  def _getEstimate(element):
    raise RuntimeError, "not implemented"




class tAnalyticSolutionErrorEstimator(tErrorEstimator):
  def __init__(self, mesh, computed_solution, analytic_solution):
    tErrorEstimator.__init__(self, mesh, computed_solution)
    self.AnalyticSolution = analytic_solution

  def _getEstimate(self, element):
    def errorFunctionL2(point, solution_func_value):
      return (self.AnalyticSolution(point) - solution_func_value) ** 2

    node_values = num.take(self.Solution, element.nodeNumbers())
    return element.getVolumeIntegralOver(errorFunctionL2, node_values)

