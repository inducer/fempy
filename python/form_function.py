import Numeric as num
import LinearAlgebra as la
import expression




# helpers ---------------------------------------------------------------------
def generateCombinations( elements, count ):
  if count == 0:
    return [ [] ]
  result = []
  for el_index in range( len( elements ) ):
    result += [ [elements[el_index]] + rest 
      for rest in generateCombinations( elements[el_index:], count - 1 ) ]
  return result

    


def multiplyUp( list_of_variables ):
  if len( list_of_variables ) == 0:
    return 1
  elif len( list_of_variables ) == 1:
    return list_of_variables[0]
  else:
    return ("*", list_of_variables[0], multiplyUp( list_of_variables[1:] ) )

def linearCombination( coefficients, expressions ):
  if len( coefficients ) == 0:
    return 0
  result = ("*", coefficients[0], expressions[0] )
  if len( coefficients ) > 1:
    result = ("+", result, linearCombination( coefficients[1:], expressions[1:] ) )
  return result





# constraints -----------------------------------------------------------------
def hasValueAt( value, point ):
  def getLinearSystemParts( expressions ):
    n = len( expressions )
    row = num.zeros( (n,), num.Float )

    var_assignments = {}
    for i in range( 0, len( point ) ):
      var_assignments[ "%d" % i ] = point[ i ]

    for i in range( 0, n ):
      row[i] = expression.evaluate( expressions[i], var_assignments )
    return row,value
  return getLinearSystemParts

def zeroAt(*point):
  return hasValueAt( 0, point )

def oneAt(*point):
  return hasValueAt( 1, point )




# form function creator -------------------------------------------------------
def makeFormFunctionExpression( order, dimensions, constraints, extra_expressions = [] ):
  variables = [ ("variable","%d" % dim ) for dim in range( 0, dimensions ) ]
  expressions = [ 1 ] + [ multiplyUp( combination )
    for ord in range( 1, order + 1 )
    for combination in generateCombinations( variables, ord ) ] + extra_expressions

  n = len( expressions )

  if n != len( constraints ):
    raise RuntimeError, "Form function has %d constraints and %d coefficients" % ( len( constraints ), n )
  
  matrix = num.zeros( (n,n), num.Float )
  rhs = num.zeros( (n,), num.Float )
  for row in range( 0, n ):
    lh_row, rh_scalar = constraints[ row ]( expressions )
    matrix[ row ] = lh_row
    rhs[ row ] = rh_scalar

  coefficients = la.solve_linear_equations( matrix, rhs )
  return expression.simplify( linearCombination( coefficients, expressions ) )
