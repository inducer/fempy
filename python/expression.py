import pattern
import string
from pattern import VAR

# This is a miniature Lisp.
# Yeah, I know. Sorry, but it is.

def evaluate( expression, variable_assignments = {} ):
  def ev( expression ):
    return pattern.switch( ruleset, expression )

  ruleset = [
    ("+",VAR,VAR), lambda x,y: ev(x)+ev(y),
    ("-",VAR,VAR), lambda x,y: ev(x)-ev(y),
    ("*",VAR,VAR), lambda x,y: ev(x)*ev(y),
    ("/",VAR,VAR), lambda x,y: ev(x)/ev(y),
    ("**",VAR,VAR), lambda x,y: ev(x)**ev(y),
    ("-",VAR) , lambda x: -ev(x),
    ("variable",VAR), lambda x: variable_assignments[x],
    VAR, lambda x: x 
  ]

  return ev( expression )




def isConstant( expression, variable_assignments = {} ):
  def isc( expression ):
    return pattern.switch( ruleset, expression )

  ruleset = [
    ("+",VAR,VAR), lambda x,y: isc(x) and isc(y),
    ("-",VAR,VAR), lambda x,y: isc(x) and isc(y),
    ("*",VAR,VAR), lambda x,y: isc(x) and isc(y),
    ("/",VAR,VAR), lambda x,y: isc(x) and isc(y),
    ("**",VAR,VAR), lambda x,y: isc(x) and isc(y),
    ("-",VAR) , lambda x: isc(x),
    ("variable",VAR), lambda x: False,
    VAR, lambda x: True
  ]

  return isc( expression )




def differentiate( expression, variable ):
  def diff( expression ):
    return pattern.switch( ruleset, expression )
  
  def diffVariable( var ):
    if var == variable:
      return 1
    else:
      return 0

  ruleset = [
    ("+",VAR,VAR), lambda x,y: ("+", diff(x), diff(y)),
    ("-",VAR,VAR), lambda x,y: ("-", diff(x), diff(y)),
    ("*",VAR,VAR), lambda x,y: ("+", ("*",diff(x),y),("*",x,diff(y))),
    ("/",VAR,VAR), lambda f,g: ("/", ("-", ("*",diff(f),g),("*",f,diff(g))), ("**", g, 2 )),
    # This assumes that the exponent is constant.
    ("**",VAR,VAR), lambda x,y: ("*",y,("*",("**",x,y-1),diff(x))),
    ("-",VAR) , lambda x: ("-",diff(x)),
    ("variable",VAR), diffVariable,
    VAR, lambda x: 0.
  ]

  return diff( expression )




def simplify( expression ):
  if isConstant( expression ):
    return evaluate( expression )

  def simplifyPlus( x, y ):
    x = simp( x )
    y = simp( y )

    if x == 0:
      return y
    if y == 0:
      return x
    return ("+",x,y)

  def simplifyMinus( x, y ):
    x = simp( x )
    y = simp( y )

    if x == 0:
      return ("-",y)
    if y == 0:
      return x
    return ("-",x,y)
  
  def simplifyTimes( x, y ):
    x = simp( x )
    y = simp( y )

    if x == 0 or y == 0:
      return 0
    if x == 1:
      return y
    if y == 1:
      return x
    return ("*",x,y)

  def simplifyDivision( x, y ):
    x = simp( x )
    y = simp( y )

    if y == 1:
      return x
    return ("/",x,y)

  def simplifyPower( x, y ):
    x = simp( x )

    if y == 0:
      return 1
    if y == 1:
      return x
    return ("**",x,y)

  def simp( expression ):
    return pattern.switch( ruleset, expression )
  
  ruleset = [
    ("+",("-",VAR),("-",VAR)), lambda x,y: ("-",("+",simplify( x ), simplify( y ))),
    ("+",VAR,("-",VAR)), lambda x,y: simplify( ("-",simplify( x ),simplify( y ))),
    ("*",-1.,VAR), lambda x: ("-",simplify( x )),
    ("*",VAR,-1.), lambda x: ("-", simplify( x )),

    ("+",VAR,VAR), simplifyPlus,
    ("-",VAR,VAR), simplifyMinus,
    ("*",VAR,VAR), simplifyTimes,
    ("/",VAR,VAR), simplifyDivision,
    # This assumes that the exponent is constant.
    ("**",VAR,VAR), simplifyPower,
    ("-",VAR) , lambda x: ("-",simp(x)),
    ("variable",VAR), lambda x: ("variable",x),
    VAR, lambda x: x
  ]

  return simp( expression )




def infixify( expression, variable_substitutions = {} ):
  determined_variables = []

  def pythonify( expr ):
    return pattern.switch( ruleset, expr )

  ruleset = [
    ("+",VAR,VAR), lambda x,y: "(%s + %s)" % ( pythonify(x), pythonify( y ) ),
    ("-",VAR,VAR), lambda x,y: "(%s - %s)" % ( pythonify(x), pythonify( y ) ),
    ("*",VAR,VAR), lambda x,y: "(%s * %s)" % ( pythonify(x), pythonify( y ) ),
    ("/",VAR,VAR), lambda x,y: "(%s / %s)" % ( pythonify(x), pythonify( y ) ),
    ("**",VAR,VAR), lambda x,y: "(%s ** %s)" % ( pythonify(x), pythonify( y ) ),
    ("-",VAR) , lambda x: "(-%s)"  % pythonify(x),
    ("variable",VAR), lambda x:"%s" % str( variable_substitutions[x] ),
    VAR, lambda x: str(x)
  ]

  return pythonify( expression )




def compile( expression, variable_substitutions = {}, variables = [] ):
  determined_variables = []

  def addVariable( var ):
    if var in variable_substitutions:
      var = variable_substitutions[var]
    if var not in determined_variables:
      determined_variables.append( var )
    return var

  def pythonify( expr ):
    return pattern.switch( ruleset, expr )

  ruleset = [
    ("+",VAR,VAR), lambda x,y: "(%s + %s)" % ( pythonify(x), pythonify( y ) ),
    ("-",VAR,VAR), lambda x,y: "(%s - %s)" % ( pythonify(x), pythonify( y ) ),
    ("*",VAR,VAR), lambda x,y: "(%s * %s)" % ( pythonify(x), pythonify( y ) ),
    ("/",VAR,VAR), lambda x,y: "(%s / %s)" % ( pythonify(x), pythonify( y ) ),
    ("**",VAR,VAR), lambda x,y: "(%s ** %s)" % ( pythonify(x), pythonify( y ) ),
    ("-",VAR) , lambda x: "(-%s)"  % pythonify(x),
    ("variable",VAR), addVariable,
    VAR, lambda x: str(x)
  ]

  pythonified = pythonify( expression )

  determined_variables.sort()
  if len( variables ) == 0 and len(determined_variables ) != 0:
    variable_str = ",".join( determined_variables )
  else:
    variable_str = ",".join( variables )
  return eval( "lambda %s:%s" % (variable_str, pythonified) )




#sample = ("+", ("variable","x"),("*", 4,("-",("variable","x")))) 

#print compile( sample )(3)
#print evaluate( sample, {"x":3} )
#print differentiate( sample, "x" )
