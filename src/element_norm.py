import tools
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
      xf_jac = element.getTransformJacobian(point)
      xf_jac_inv = la.inverse(xf_jac)
      grad = element.getFormFunctionCombinationGradient(node_values, point)
      twisted_grad = num.matrixmultiply(num.transpose(xf_jac_inv), grad)

      return ((analytic_solution(real_point) - solution_func_value) ** 2 \
        + tools.norm2squared(twisted_grad - grad_analytic_solution(real_point)))

    node_values = num.take(solution_vector, element.nodeNumbers())

    return element.getVolumeIntegralOver(errorFunctionH1, node_values)
  return tools.tFunctionValueCache(result_func)

def makeEnergyErrorNormSquared(grad_analytic_solution, solution_vector):
  def result_func(element):
    def errorFunctionH1(point, solution_func_value):
      real_point = element.transformToReal(point)
      xf_jac = element.getTransformJacobian(point)
      xf_jac_inv = la.inverse(xf_jac)
      grad = element.getFormFunctionCombinationGradient(node_values, point)
      twisted_grad = num.matrixmultiply(num.transpose(xf_jac_inv), grad)

      return tools.norm2squared(twisted_grad - grad_analytic_solution(real_point))

    node_values = num.take(solution_vector, element.nodeNumbers())
    approx_grad_solution = element.getFormFunctionCombinationGradient(node_values)

    return element.getVolumeIntegralOver(errorFunctionH1, node_values)
  return tools.tFunctionValueCache(result_func)
