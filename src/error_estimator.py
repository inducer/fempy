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

