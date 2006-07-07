import pylinear.array as num
import pylinear.linear_algebra as la
import expression_operators as eo
import expression




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
def hasValueAt(value, point):
    def get_linear_system_parts(expressions):
        n = len(expressions)
        row = num.zeros((n,), num.Float)
        
        var_assignments = {}
        for i in range(0, len(point)):
            var_assignments[ "%d" % i ] = point[ i ]

        for i in range(0, n):
            row[i] = expression.evaluate(expressions[i], var_assignments)
        return row,value
    return get_linear_system_parts

def zero_at(*point):
  return has_value_at(0, point)

def one_at(*point):
  return has_value_at(1, point)




# form function creator -------------------------------------------------------
def make_form_function_expression(order, dimensions, constraints, extra_expressions=[]):
    variables = [ (eo.VARIABLE,"%d" % dim) for dim in range(0, dimensions) ]
    expressions = [ 1 ] + [ 
        expression.multiplyUp(combination)
        for current_ord in range(1, order + 1)
        for combination in generateCombinations(variables, current_ord) ] \
        + extra_expressions

    n = len(expressions)

    if n != len(constraints):
        raise RuntimeError, "Form function has %d constraints and %d coefficients" % (len(constraints), n)
  
    matrix = num.zeros((n,n), num.Float)
    rhs = num.zeros((n,), num.Float)
    for row in range(0, n):
        lh_row, rh_scalar = constraints[row](expressions)
        matrix[row] = lh_row
        rhs[row] = rh_scalar

    coefficients = la.solve_linear_equations(matrix, rhs)
    return expression.simplify(expression.linearCombination(coefficients, expressions))




def makeFormFunctions(order, points, extra_expressions=[]):
    all_constraints = []
    for i in range(len(points)):
        all_constraints.append([zeroAt(*point) for point in points])
        all_constraints[i][i] = oneAt(*points[i])
    return [makeFormFunctionExpression(
        order, len(points[0]), constraint, 
        extra_expressions) for constraint in all_constraints]

