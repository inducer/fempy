"""This is heavily inspired by OCaml's pattern matching. While its Python equivalent
is fairly slow, it is still useful."""

import types

class VAR:
  pass

def match(pattern, value):
  if pattern == VAR:
    return [value]
  if type(pattern) != type(value):
    return None
  if pattern == value:
    return []
  if type(pattern) in [ types.TupleType, types.ListType ] :
    if len(pattern) != len(value):
      return None
    result_vars = []
    for index in range(0, len(pattern)):
      result = match(pattern[ index ], value[ index ])
      if result is None:
	return None
      else:
	result_vars.extend(result)
    return result_vars




def switch(ruleset, value):
  for index in range(0, len(ruleset), 2):
    variables = match(ruleset[ index ], value)
    if variables is not None:
      return ruleset[ index + 1 ](*variables)
  raise RuntimeError, "Value did not match a pattern"
