import pyanglemetal
from pyanglemetal import tTriangulationParameters




def triangulate( in_parameters, out_parameters, verbose = False,
    refinement_func = None):
  opts = "pzq"
  if verbose:
    opts += "V"
  if refinement_func is not None:
    opts += "u"
  pyanglemetal.triangulate( 
      opts, in_parameters, out_parameters, 
      tTriangulationParameters(),
      refinement_func )

