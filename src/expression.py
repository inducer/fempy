import pattern
from pattern import VAR
import types
import pylinear.matrices as num
import math




def ispos(x):
  value = ev(x)
  if value > 0:
    return 1
  else:
    return 0




def matchRuleSet(ruleset, expression):
  try:
    return ruleset[expression[0], len(expression)-1](*expression[1:])
  except TypeError:
    return ruleset[None](expression)




def evaluate(expression, variable_assignments = {}):
  def ev(expression):
    return matchRuleSet(ruleset, expression)

  ruleset = {
    ("+",2): lambda x,y: ev(x)+ev(y),
    ("-",2): lambda x,y: ev(x)-ev(y),
    ("*",2): lambda x,y: ev(x)*ev(y),
    ("/",2): lambda x,y: ev(x)/ev(y),
    ("**",2): lambda x,y: ev(x)**ev(y),
    ("-",1): lambda x: -ev(x),
    ("sin",1): lambda x: math.sin(ev(x)),
    ("cos",1): lambda x: math.cos(ev(x)),
    ("tan",1): lambda x: math.tan(ev(x)),
    ("ispos",1): lambda x: ispos(ev(x)),
    ("variable",1): lambda x: variable_assignments[x],
    None: lambda x: x 
  }

  return ev(expression)




def substitute(expression, variable_assignments = {}):
  def subs(expression):
    return pattern.switch(ruleset, expression)

  def substituteVariable(var_id):
    try:
      return variable_assignments[var_id]
    except KeyError:
      return ("variable", var_id)

  ruleset = [
    ("+",VAR,VAR), lambda x,y: ("+",subs(x),subs(y)),
    ("-",VAR,VAR), lambda x,y: ("-",subs(x),subs(y)),
    ("*",VAR,VAR), lambda x,y: ("*",subs(x),subs(y)),
    ("/",VAR,VAR), lambda x,y: ("/",subs(x),subs(y)),
    ("**",VAR,VAR), lambda x,y: ("**",subs(x),subs(y)),
    ("-",VAR) , lambda x: ("-",subs(x)),
    ("sin",VAR), lambda x: ("sin",subs(x)),
    ("cos",VAR), lambda x: ("cos",subs(x)),
    ("tan",VAR), lambda x: ("tan",subs(x)),
    ("ispos",VAR), lambda x: ("ispos",subs(x)),
    ("variable",VAR), substituteVariable,
    VAR, lambda x: x 
  ]

  return subs(expression)




def isConstant(expression, wrt = None):
  """Returns whether the expression is constant, i.e. does not 
  depend on variables.

  If wrt (short for with respect to) is None, the occurrence of any 
  variable makes the expression non-constant.

  If wrt is a list of variable identifiers, only those variables 
  are assumed "non-constant".
  """
  def isc(expression):
    return pattern.switch(ruleset, expression)

  def isVarConstant(x):
    if wrt:
      return x not in wrt
    else:
      return False

  ruleset = [
    ("+",VAR,VAR), lambda x,y: isc(x) and isc(y),
    ("-",VAR,VAR), lambda x,y: isc(x) and isc(y),
    ("*",VAR,VAR), lambda x,y: isc(x) and isc(y),
    ("/",VAR,VAR), lambda x,y: isc(x) and isc(y),
    ("**",VAR,VAR), lambda x,y: isc(x) and isc(y),
    ("-",VAR) , lambda x: isc(x),
    ("sin",VAR), lambda x: isc(x),
    ("cos",VAR), lambda x: isc(x),
    ("tan",VAR), lambda x: isc(x),
    ("ispos",VAR), lambda x: isc(x),
    ("variable",VAR), isVarConstant,
    VAR, lambda x: True
  ]

  return isc(expression)




def differentiate(expression, variable):
  def diff(expression):
    if isConstant(expression, [variable]):
      return 0
    return matchRuleSet(ruleset, expression)
  
  def diffVariable(var):
    if var == variable:
      return 1
    else:
      return 0

  def ispos_diff(x):
    raise RuntimeError, "ispos is not differentiable"

  def diffMultiplication(f,g):
    if isConstant(f, [variable]):
      return ("*",f,diff(g))
    elif isConstant(g, [variable]):
      return ("*",g,diff(f))
    else:
      return ("+", ("*",diff(f),g),("*",f,diff(g)))

  def diffDivision(f,g):
    if isConstant(f, [variable]):
      return ("/", ("-",("*",f,diff(g))), ("**", g, 2))
    elif isConstant(g, [variable]):
      return ("/",diff(f),g)
    else:
      return ("/", ("-", ("*",diff(f),g),("*",f,diff(g))), ("**", g, 2))

  ruleset = {
    ("+",2): lambda x,y: ("+", diff(x), diff(y)),
    ("-",2): lambda x,y: ("-", diff(x), diff(y)),
    ("*",2): diffMultiplication,
    ("/",2): diffDivision,
    # This assumes that the exponent is constant.
    ("**",2): lambda x,y: ("*",y,("*",("**",x,y-1),diff(x))),
    ("-",1): lambda x: ("-",diff(x)),
    ("sin",1): lambda x: ("*", diff(x), ("cos",x)),
    ("cos",1): lambda x: ("*", ("-", diff(x)), ("sin",x)),
    ("tan",1): lambda x: ("*", diff(x), ("+", 1, ("**", ("tan", x), 2))),
    ("ispos",1): ispos_diff,
    ("variable",1): diffVariable,
    None: lambda x: 0.
  }

  return diff(expression)




def simplify(expression):
  if isConstant(expression):
    return evaluate(expression)

  def simplifyPlus(x, y):
    x = simp(x)
    y = simp(y)

    if x == 0:
      return y
    if y == 0:
      return x
    return ("+",x,y)

  def simplifyMinus(x, y):
    x = simp(x)
    y = simp(y)

    if x == 0:
      return ("-",y)
    if y == 0:
      return x
    return ("-",x,y)
  
  def simplifyTimes(x, y):
    x = simp(x)
    y = simp(y)

    if x == 0 or y == 0:
      return 0
    if x == 1:
      return y
    if y == 1:
      return x
    return ("*",x,y)

  def simplifyDivision(x, y):
    x = simp(x)
    y = simp(y)

    if y == 1:
      return x
    return ("/",x,y)

  def simplifyPower(x, y):
    x = simp(x)

    if y == 0:
      return 1
    if y == 1:
      return x
    return ("**",x,y)

  def simp(expression):
    return matchRuleSet(ruleset, expression)
  
  ruleset = {
    #("+",("-",VAR),("-",VAR)), lambda x,y: ("-",("+",simplify(x), simplify(y))),
    #("+",VAR,("-",VAR)), lambda x,y: simplify(("-",simplify(x),simplify(y))),
    #("*",-1.,VAR), lambda x: ("-",simplify(x)),
    #("*",VAR,-1.), lambda x: ("-", simplify(x)),

    ("+",2): simplifyPlus,
    ("-",2): simplifyMinus,
    ("*",2): simplifyTimes,
    ("/",2): simplifyDivision,
    # This assumes that the exponent is constant.
    ("**",2): simplifyPower,
    ("-",1): lambda x: ("-",simp(x)),
    ("sin",1): lambda x: ("sin",simp(x)),
    ("cos",1): lambda x: ("cos",simp(x)),
    ("tan",1): lambda x: ("tan",simp(x)),
    ("ispos",1): lambda x: ("ispos",simp(x)),
    ("variable",1): lambda x: ("variable",x),
    None: lambda x: x
  }

  return simp(expression)




def infixify(expression, variable_substitutions = {}):
  determined_variables = []

  def substitute(var):
    try:
      return variable_substitutions[var]
    except KeyError:
      return "V"+str(var)

  def pythonify(expr):
    return pattern.switch(ruleset, expr)

  ruleset = [
    ("+",VAR,VAR), lambda x,y: "(%s + %s)" % (pythonify(x), pythonify(y)),
    ("-",VAR,VAR), lambda x,y: "(%s - %s)" % (pythonify(x), pythonify(y)),
    ("*",VAR,VAR), lambda x,y: "(%s * %s)" % (pythonify(x), pythonify(y)),
    ("/",VAR,VAR), lambda x,y: "(%s / %s)" % (pythonify(x), pythonify(y)),
    ("**",VAR,VAR), lambda x,y: "(%s ** %s)" % (pythonify(x), pythonify(y)),
    ("-",VAR) , lambda x: "(-%s)"  % pythonify(x),
    ("sin",VAR) , lambda x: "sin(%s)"  % pythonify(x),
    ("cos",VAR) , lambda x: "cos(%s)"  % pythonify(x),
    ("tan",VAR) , lambda x: "tan(%s)"  % pythonify(x),
    ("ispos",VAR) , lambda x: "ispos(%s)"  % pythonify(x),
    ("variable",VAR), lambda x:"%s" % str(substitute(x)),
    VAR, lambda x: str(x)
  ]

  return pythonify(expression)




def compile(expression, variable_substitutions = {}, variables = []):
  determined_variables = []

  def addVariable(var):
    if var in variable_substitutions:
      var = variable_substitutions[var]
    if var not in determined_variables:
      determined_variables.append(var)
    return var

  def pythonify(expr):
    return pattern.switch(ruleset, expr)

  ruleset = [
    ("+",VAR,VAR), lambda x,y: "(%s + %s)" % (pythonify(x), pythonify(y)),
    ("-",VAR,VAR), lambda x,y: "(%s - %s)" % (pythonify(x), pythonify(y)),
    ("*",VAR,VAR), lambda x,y: "(%s * %s)" % (pythonify(x), pythonify(y)),
    ("/",VAR,VAR), lambda x,y: "(%s / %s)" % (pythonify(x), pythonify(y)),
    ("**",VAR,VAR), lambda x,y: "(%s ** %s)" % (pythonify(x), pythonify(y)),
    ("-",VAR) , lambda x: "(-%s)"  % pythonify(x),
    ("sin",VAR) , lambda x: "math.sin(%s)"  % pythonify(x),
    ("cos",VAR) , lambda x: "math.cos(%s)"  % pythonify(x),
    ("tan",VAR) , lambda x: "math.tan(%s)"  % pythonify(x),
    ("ispos",VAR) , lambda x: "ispos(%s)"  % pythonify(x),
    ("variable",VAR), addVariable,
    VAR, lambda x: str(x)
  ]

  pythonified = pythonify(expression)

  determined_variables.sort()
  if len(variables) == 0 and len(determined_variables) != 0:
    variable_str = ",".join(determined_variables)
  else:
    variable_str = ",".join(variables)
  return eval("lambda %s:%s" % (variable_str, pythonified))




# tools -----------------------------------------------------------------------
def compileScalarField(expr, dimension = 2):
  """Compiles a function that has variables named "0", "1", "2" up to
  dimension into a function that has one sequence as its input.
  """
  substitutions = {}
  for i in range(dimension):
    substitutions[ "%d" % i ] = "point[%d]" % i

  return compile(expr, substitutions, [ "point" ])




def assembleVectorFunction(function_list):
  return lambda x: num.array([ func(x) for func in function_list ])




def assembleMatrixFunction(function_list):
  return lambda x: num.array([ [ func(x) for func in flist ] for flist in function_list ])




def compileVectorField(expr_list, preimage_dim = 2):
  return assembleVectorFunction([compileScalarField(expr, preimage_dim) for expr in expr_list])




def compileMatrixFunction(expr_list, preimage_dim = 2):
  return assembleMatrixFunction(
    [
    [
    compileScalarField(expr, preimage_dim)
    for expr in outer]
    for outer in expr_list])





def sumUp(list_of_expressions):
  if len(list_of_expressions) == 0:
    return 0
  elif len(list_of_expressions) == 1:
    return list_of_expressions[0]
  else:
    return ("+", list_of_expressions[0], sumUp(list_of_expressions[1:]))




def multiplyUp(list_of_expressions):
  if len(list_of_expressions) == 0:
    return 1
  elif len(list_of_expressions) == 1:
    return list_of_expressions[0]
  else:
    return ("*", list_of_expressions[0], multiplyUp(list_of_expressions[1:]))




def linearCombination(coefficients, expressions):
  if len(coefficients) == 0:
    return 0
  result = ("*", coefficients[0], expressions[0])
  if len(coefficients) > 1:
    result = ("+", result, linearCombination(coefficients[1:], expressions[1:]))
  return result




def grad(expression, variables):
  return [simplify(differentiate(expression, var)) for var in variables]




def jacobian(expression_list, variables):
  return [grad(expr, variables) for expr in expression_list]




def laplace(expression, variables):
  return sumUp([differentiate(differentiate(expression,var),var) for var in variables])





