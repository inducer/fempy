import tools
import mesh_function
import pylinear.matrices as num
import pylinear.linear_algebra as la




def makeL2ErrorNormSquared(analytic_solution, solution_vector):
  def result_func(element):
    def errorFunctionL2(point, solution_func_value):
      return (analytic_solution(element.transformToReal(point)) - solution_func_value) ** 2

    node_values = num.take(solution_vector, element.nodeNumbers())
    return element.getVolumeIntegralOver(errorFunctionL2, node_values)
  return tools.tFunctionValueCache(result_func)

def makeH1ErrorNormSquared(analytic_solution, grad_analytic_solution, solution_vector):
  def result_func(element):
    def errorFunctionH1(point, solution_func_value):
      real_point = element.transformToReal(point)
      grad = mesh_function.getRealGradientOnElement(element, solution_vector, point)

      return ((analytic_solution(real_point) - solution_func_value) ** 2 \
        + tools.norm2squared(grad - grad_analytic_solution(real_point)))

    node_values = num.take(solution_vector, element.nodeNumbers())
    return element.getVolumeIntegralOver(errorFunctionH1, node_values)

  return tools.tFunctionValueCache(result_func)

def makeEnergyErrorNormSquared(grad_analytic_solution, solution_vector):
  def result_func(element):
    def errorFunctionH1(point, solution_func_value):
      real_point = element.transformToReal(point)
      grad = mesh_function.getRealGradientOnElement(element, solution_vector, point)

      return tools.norm2squared(grad - grad_analytic_solution(real_point))

    node_values = num.take(solution_vector, element.nodeNumbers())
    return element.getVolumeIntegralOver(errorFunctionH1, node_values)

  return tools.tFunctionValueCache(result_func)
