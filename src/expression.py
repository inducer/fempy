import pattern
import pylinear.matrices as num
import math
import sets
import expression_operators as eo

VAR = pattern.VAR




def ispos(x):
    if x > 0:
        return 1
    else:
        return 0




def _matchRuleSet(ruleset, expression):
    try:
        return ruleset[expression[0], len(expression)-1](*expression[1:])
    except TypeError:
        return ruleset[None](expression)




def evaluate(expression, variable_assignments = {}):
    def ev(expression):
        return _matchRuleSet(ruleset, expression)

    ruleset = {
        (eo.PLUS,2): lambda x,y: ev(x)+ev(y),
        (eo.MINUS,2): lambda x,y: ev(x)-ev(y),
        (eo.TIMES,2): lambda x,y: ev(x)*ev(y),
        (eo.DIVIDE,2): lambda x,y: ev(x)/ev(y),
        (eo.POWER,2): lambda x,y: ev(x)**ev(y),
        (eo.NEG,1): lambda x: -ev(x),
        (eo.SIN,1): lambda x: math.sin(ev(x)),
        (eo.COS,1): lambda x: math.cos(ev(x)),
        (eo.TAN,1): lambda x: math.tan(ev(x)),
        (eo.LOG,1): lambda x: math.log(ev(x)),
        (eo.EXP,1): lambda x: math.exp(ev(x)),
        (eo.ISPOS,1): lambda x: ispos(ev(x)),
        (eo.VARIABLE,1): lambda x: variable_assignments[x],
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
            return (eo.VARIABLE, var_id)

    ruleset = [
        ("+",VAR,VAR), lambda x,y: (eo.PLUS,subs(x),subs(y)),
        ("-",VAR,VAR), lambda x,y: (eo.MINUS,subs(x),subs(y)),
        ("*",VAR,VAR), lambda x,y: (eo.TIMES,subs(x),subs(y)),
        ("/",VAR,VAR), lambda x,y: (eo.DIVIDE,subs(x),subs(y)),
        ("**",VAR,VAR), lambda x,y: (eo.POWER,subs(x),subs(y)),
        ("negate",VAR) , lambda x: (eo.NEG,subs(x)),
        ("sin",VAR), lambda x: (eo.SIN,subs(x)),
        ("cos",VAR), lambda x: (eo.COS,subs(x)),
        ("tan",VAR), lambda x: (eo.TAN,subs(x)),
        ("log",VAR), lambda x: (eo.LOG,subs(x)),
        ("exp",VAR), lambda x: (eo.EXP,subs(x)),
        ("ispos",VAR), lambda x: (eo.ISPOS,subs(x)),
        ("variable",VAR), substituteVariable,
        VAR, lambda x: x 
        ]

    substituted = subs(expression)
    if substituted != expression:
        return substituted
    else:
        # the substituted copy will be gc'ed.
        return expression




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
        (eo.NEG,VAR) , lambda x: isc(x),
        ("sin",VAR), lambda x: isc(x),
        ("cos",VAR), lambda x: isc(x),
        ("tan",VAR), lambda x: isc(x),
        ("log",VAR), lambda x: isc(x),
        ("exp",VAR), lambda x: isc(x),
        ("ispos",VAR), lambda x: isc(x),
        ("variable",VAR), isVarConstant,
        VAR, lambda x: True
        ]

    return isc(expression)




def differentiate(expression, variable):
    def diff(expression):
        if isConstant(expression, [variable]):
            return 0
        return _matchRuleSet(ruleset, expression)
  
    def diffVariable(var):
        if var == variable:
            return 1
        else:
            return 0

    def ispos_diff(x):
        raise RuntimeError, "ispos is not differentiable"

    def diffMultiplication(f,g):
        if isConstant(f, [variable]):
            return (eo.TIMES,f,diff(g))
        elif isConstant(g, [variable]):
            return (eo.TIMES,g,diff(f))
        else:
            return (eo.PLUS, (eo.TIMES,diff(f),g),(eo.TIMES,f,diff(g)))

    def diffDivision(f,g):
        if isConstant(f, [variable]):
            return (eo.DIVIDE, (eo.NEG,(eo.TIMES,f,diff(g))), 
                    (eo.POWER, g, 2))
        elif isConstant(g, [variable]):
            return (eo.DIVIDE,diff(f),g)
        else:
            return (eo.DIVIDE, (eo.MINUS, (eo.TIMES,diff(f),g),(eo.TIMES,f,diff(g))), 
                    (eo.POWER, g, 2))

    ruleset = {
        (eo.PLUS,2): lambda x,y: (eo.PLUS, diff(x), diff(y)),
        (eo.MINUS,2): lambda x,y: (eo.MINUS, diff(x), diff(y)),
        (eo.TIMES,2): diffMultiplication,
        (eo.DIVIDE,2): diffDivision,
        # This assumes that the exponent is constant.
        (eo.POWER,2): lambda x,y: (eo.TIMES,y,(eo.TIMES,(eo.POWER,x,y-1),diff(x))),
        (eo.NEG,1): lambda x: (eo.NEG,diff(x)),
        (eo.SIN,1): lambda x: (eo.TIMES, diff(x), (eo.COS,x)),
        (eo.COS,1): lambda x: (eo.TIMES, (eo.NEG, diff(x)), (eo.SIN,x)),
        (eo.TAN,1): lambda x: (eo.TIMES, diff(x), (eo.PLUS, 1, (eo.POWER, (eo.TAN, x), 2))),
        (eo.LOG,1): lambda x: (eo.TIMES, diff(x), (eo.DIVIDE, 1, x)),
        (eo.EXP,1): lambda x: (eo.TIMES, diff(x), (eo.EXP, x)),
        (eo.ISPOS,1): ispos_diff,
        (eo.VARIABLE,1): diffVariable,
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
        return (eo.PLUS,x,y)

    def simplifyMinus(x, y):
        x = simp(x)
        y = simp(y)
        
        if x == 0:
            return (eo.NEG,y)
        if y == 0:
            return x
        return (eo.MINUS,x,y)
  
    def simplifyTimes(x, y):
        x = simp(x)
        y = simp(y)

        if x == 0 or y == 0:
            return 0
        if x == 1:
            return y
        if y == 1:
            return x
        return (eo.TIMES,x,y)

    def simplifyDivision(x, y):
        x = simp(x)
        y = simp(y)
        
        if y == 1:
            return x
        return (eo.DIVIDE,x,y)

    def simplifyPower(x, y):
        x = simp(x)
        
        if y == 0:
            return 1
        if y == 1:
            return x
        return (eo.POWER,x,y)

    def simp(expression):
        return _matchRuleSet(ruleset, expression)
  
    ruleset = {
        #("+",("-",VAR),("-",VAR)), lambda x,y: ("-",("+",simplify(x), simplify(y))),
        #("+",VAR,("-",VAR)), lambda x,y: simplify(("-",simplify(x),simplify(y))),
        #("*",-1.,VAR), lambda x: ("-",simplify(x)),
        #("*",VAR,-1.), lambda x: ("-", simplify(x)),
        
        (eo.PLUS,2): simplifyPlus,
        (eo.MINUS,2): simplifyMinus,
        (eo.TIMES,2): simplifyTimes,
        (eo.DIVIDE,2): simplifyDivision,
        # This assumes that the exponent is constant.
        (eo.POWER,2): simplifyPower,
        (eo.NEG,1): lambda x: (eo.NEG,simp(x)),
        (eo.SIN,1): lambda x: (eo.SIN,simp(x)),
        (eo.COS,1): lambda x: (eo.COS,simp(x)),
        (eo.TAN,1): lambda x: (eo.TAN,simp(x)),
        (eo.LOG,1): lambda x: (eo.LOG,simp(x)),
        (eo.EXP,1): lambda x: (eo.EXP,simp(x)),
        (eo.ISPOS,1): lambda x: (eo.ISPOS,simp(x)),
        (eo.VARIABLE,1): lambda x: (eo.VARIABLE,x),
        None: lambda x: x
        }
    
    return simp(expression)




def infixify(expression, variable_substitutions = {}):
    def substitute(var):
        try:
            return variable_substitutions[var]
        except KeyError:
            return "V"+str(var)

    def pythonify(expr):
        return _matchRuleSet(ruleset, expr)

    ruleset = {
        (eo.PLUS,2): lambda x,y: "(%s + %s)" % (pythonify(x), pythonify(y)),
        (eo.MINUS,2): lambda x,y: "(%s - %s)" % (pythonify(x), pythonify(y)),
        (eo.TIMES,2): lambda x,y: "(%s * %s)" % (pythonify(x), pythonify(y)),
        (eo.DIVIDE,2): lambda x,y: "(%s / %s)" % (pythonify(x), pythonify(y)),
        (eo.POWER,2): lambda x,y: "(%s ** %s)" % (pythonify(x), pythonify(y)),
        (eo.NEG,1): lambda x: "(-%s)"  % pythonify(x),
        (eo.SIN,1): lambda x: "sin(%s)"  % pythonify(x),
        (eo.COS,1): lambda x: "cos(%s)"  % pythonify(x),
        (eo.TAN,1): lambda x: "tan(%s)"  % pythonify(x),
        (eo.LOG,1): lambda x: "log(%s)"  % pythonify(x),
        (eo.EXP,1): lambda x: "exp(%s)"  % pythonify(x),
        (eo.ISPOS,1): lambda x: "ispos(%s)"  % pythonify(x),
        (eo.VARIABLE,1): lambda x:"%s" % str(substitute(x)),
        None: lambda x: str(x)
        }    

    return pythonify(expression)




class tCompiledExpression:
    def __init__(self, expression, variable_substitutions = {}, variables = []):
        self.Expression = expression
        self.VariableSubstitutions = variable_substitutions.copy()
        self.Variables = variables[:]
        self.__compile__()

    def __compile__(self):
        used_variables = sets.Set()

        def addVariable(var):
            try:
                var = self.VariableSubstitutions[var]
            except:
                pass
            used_variables.add(var)
            return var

        def pythonify(expr):
            return _matchRuleSet(ruleset, expr)

        ruleset = {
            (eo.PLUS,2): lambda x,y: "(%s + %s)" % (pythonify(x), pythonify(y)),
            (eo.MINUS,2): lambda x,y: "(%s - %s)" % (pythonify(x), pythonify(y)),
            (eo.TIMES,2): lambda x,y: "(%s * %s)" % (pythonify(x), pythonify(y)),
            (eo.DIVIDE,2): lambda x,y: "(%s / %s)" % (pythonify(x), pythonify(y)),
            (eo.POWER,2): lambda x,y: "(%s ** %s)" % (pythonify(x), pythonify(y)),
            (eo.NEG,1): lambda x: "(-%s)"  % pythonify(x),
            (eo.SIN,1): lambda x: "math.sin(%s)"  % pythonify(x),
            (eo.COS,1): lambda x: "math.cos(%s)"  % pythonify(x),
            (eo.TAN,1): lambda x: "math.tan(%s)"  % pythonify(x),
            (eo.LOG,1): lambda x: "math.log(%s)"  % pythonify(x),
            (eo.EXP,1): lambda x: "math.exp(%s)"  % pythonify(x),
            (eo.ISPOS,1): lambda x: "ispos(%s)"  % pythonify(x),
            (eo.VARIABLE,1): addVariable,
            None: lambda x: str(x)
            }

        pythonified = pythonify(self.Expression)

        used_variables = list(used_variables)
        used_variables.sort()
        if len(self.Variables) == 0 and len(used_variables) != 0:
            variable_str = ",".join(used_variables)
        else:
            variable_str = ",".join(self.Variables)
        self.__call__ = eval("lambda %s:%s" % (variable_str, pythonified))
    
    def __getinitargs__(self):
        return self.Expression, self.VariableSubstitutions, self.Variables
        
    def __getstate__(self):
        return None

    def __setstate__(self, state):
        pass

compile = tCompiledExpression


    

# tools -----------------------------------------------------------------------
def compileScalarField(expr, dimension = 2):
    """Compiles a function that has variables named "0", "1", "2" up to
    dimension into a function that has one sequence as its input.
    """
    substitutions = {}
    for i in range(dimension):
        substitutions[ "%d" % i ] = "point[%d]" % i
        
    return compile(expr, substitutions, [ "point" ])




class tVectorFunction:
    def __init__(self, function_list):
        self.FunctionList = function_list

    def __call__(self, x):
        return num.array([ func(x) for func in self.FunctionList ])

def assembleVectorFunction(function_list):
    return tVectorFunction(function_list)




class tMatrixFunction:
    def __init__(self, function_list):
        self.FunctionList = function_list

    def __call__(self, x):
        return num.array([ [ func(x) for func in flist ] for flist in self.FunctionList ])

def assembleMatrixFunction(function_list):
    return tMatrixFunction(function_list)




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
        return (eo.PLUS, list_of_expressions[0], sumUp(list_of_expressions[1:]))




def multiplyUp(list_of_expressions):
    if len(list_of_expressions) == 0:
        return 1
    elif len(list_of_expressions) == 1:
        return list_of_expressions[0]
    else:
        return (eo.TIMES, list_of_expressions[0], multiplyUp(list_of_expressions[1:]))




def linearCombination(coefficients, expressions):
    if len(coefficients) == 0:
        return 0
    result = (eo.TIMES, coefficients[0], expressions[0])
    if len(coefficients) > 1:
        result = (eo.PLUS, result, linearCombination(coefficients[1:], expressions[1:]))
    return result




def grad(expression, variables):
    return [simplify(differentiate(expression, var)) for var in variables]




def jacobian(expression_list, variables):
    return [grad(expr, variables) for expr in expression_list]




def laplace(expression, variables):
    return sumUp([differentiate(differentiate(expression,var),var) for var in variables])

