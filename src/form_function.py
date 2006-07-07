import pylinear.array as num
import pylinear.operation as op
import pymbolic




# helpers ---------------------------------------------------------------------
def generate_combinations(elements, count):
    if count == 0:
        return [ [] ]
    result = []
    for el_index in range(len(elements)):
        result += [[elements[el_index]] + rest 
                   for rest in generate_combinations(elements[el_index:], count - 1) ]
    return result

    


# constraints -----------------------------------------------------------------
def has_value_at(value, point):
    def get_linear_system_parts(expressions, variable):
        n = len(expressions)
        row = num.zeros((n,), num.Float)
        
        for i in range(0, n):
            row[i] = pymbolic.evaluate(expressions[i], {variable.name: point})
        return row, value
    return get_linear_system_parts

def zero_at(*point):
  return has_value_at(0, point)

def one_at(*point):
  return has_value_at(1, point)




# form function creator -------------------------------------------------------
def make_form_function_expression(order, dimensions, constraints, extra_expressions=[],
                                  variable=pymbolic.var("x")):
    variables = [pymbolic.subscript(variable, pymbolic.const(i))
                 for i in range(dimensions)]
                                    
    expressions = [ 
        pymbolic.product(*combination)
        for current_ord in range(0, order + 1)
        for combination in generate_combinations(variables, current_ord) ] \
        + extra_expressions

    n = len(expressions)

    if n != len(constraints):
        raise RuntimeError, "Form function has %d constraints and %d coefficients" % (len(constraints), n)
  
    matrix = num.zeros((n,n), num.Float)
    rhs = num.zeros((n,), num.Float)
    for row in range(0, n):
        lh_row, rh_scalar = constraints[row](expressions, variable)
        matrix[row] = lh_row
        rhs[row] = rh_scalar

    coefficients = matrix <<num.solve>> rhs
    return pymbolic.simplify(pymbolic.linear_combination(coefficients, expressions))




def make_form_functions(order, points, extra_expressions=[], 
                        variable=pymbolic.var("x")):
    all_constraints = []
    for i in range(len(points)):
        all_constraints.append([zero_at(*point) for point in points])
        all_constraints[i][i] = one_at(*points[i])
    return [make_form_function_expression(
        order, len(points[0]), constraint, 
        extra_expressions, variable) for constraint in all_constraints]

