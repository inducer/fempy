import pylinear.array as num
import pylinear.operation as op

import tools




def makeL2ErrorNormSquared(analytic_solution, mesh_func):
  def result_func(el):
    def errorFunctionL2(point, solution_func_value):
      return abs(analytic_solution(el.transformToReal(point)) - solution_func_value) ** 2

    node_values = num.array([mesh_func[node] for node in el.nodes()])
    return el.getVolumeIntegralOver(errorFunctionL2, node_values)
  return tools.FunctionValueCache(result_func)

def makeH1ErrorNormSquared(analytic_solution, grad_analytic_solution, mesh_func):
  def result_func(el):
    def errorFunctionH1(point, solution_func_value):
      real_point = element.transformToReal(point)
      grad = mesh_func.getRealGradientOnElement(el, point)

      return (abs(analytic_solution(real_point) - solution_func_value) ** 2 \
        + op.norm_2_squared(grad - grad_analytic_solution(real_point)))

    node_values = num.array([mesh_func[node] for node in el.nodes()])
    return el.getVolumeIntegralOver(errorFunctionH1, node_values)

  return tools.FunctionValueCache(result_func)

def makeEnergyErrorNormSquared(grad_analytic_solution, mesh_func):
  def result_func(el):
    def errorFunctionH1(point, solution_func_value):
      real_point = el.transformToReal(point)
      grad = mesh_func.getRealGradientOnElement(el, point)

      return op.norm_2_squared(grad - grad_analytic_solution(real_point))

    node_values = num.array([mesh_func[node] for node in el.nodes()])
    return el.getVolumeIntegralOver(errorFunctionH1, node_values)

  return tools.FunctionValueCache(result_func)
