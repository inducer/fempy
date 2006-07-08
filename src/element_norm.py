import pylinear.array as num
import pylinear.computation as comp

import pytools




def make_l2_error_norm_squared(analytic_solution, mesh_func):
    def result_func(el):
        def error_function_l2(point, solution_func_value):
            return abs(analytic_solution(el.transform_to_real(point)) - solution_func_value) ** 2

        node_values = num.array([mesh_func[node] for node in el.nodes()])
        return el.get_volume_integral_over(error_function_l2, node_values)
    return pytools.FunctionValueCache(result_func)

def make_h1_error_norm_squared(analytic_solution, grad_analytic_solution, mesh_func):
    def result_func(el):
        def error_function_h1(point, solution_func_value):
            real_point = element.transform_to_real(point)
            grad = mesh_func.get_real_gradient_on_element(el, point)

            return (abs(analytic_solution(real_point) - solution_func_value) ** 2 \
                    + comp.norm_2_squared(grad - grad_analytic_solution(real_point)))

        node_values = num.array([mesh_func[node] for node in el.nodes()])
        return el.get_volume_integral_over(error_function_h1, node_values)

    return pytools.FunctionValueCache(result_func)

def make_energy_error_norm_squared(grad_analytic_solution, mesh_func):
    def result_func(el):
        def error_function_h1(point, solution_func_value):
            real_point = el.transform_to_real(point)
            grad = mesh_func.get_real_gradient_on_element(el, point)

            return comp.norm_2_squared(grad - grad_analytic_solution(real_point))

        node_values = num.array([mesh_func[node] for node in el.nodes()])
        return el.get_volume_integral_over(error_function_h1, node_values)

    return pytools.FunctionValueCache(result_func)
