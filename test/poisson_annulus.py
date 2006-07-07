import math

import test_tools

import pymbolic
import pymbolic.functions

import fempy
import fempy.mesh
import fempy.geometry

outer_radius = 5
inner_radius = 0.1
outer_value = 0
inner_value = 10

c1 = (outer_value - inner_value)/(math.log(outer_radius/inner_radius))
c2 = inner_value - c1 * math.log(inner_radius)

expr_r = pymbolic.parse("(x[0]**2+x[1]**2)**0.5")
sol = c1 * pymbolic.functions.log(expr_r) + c2

mesh = fempy.mesh.TwoDimensionalMesh(
  fempy.geometry.get_annulus_geometry(5, 0.1, use_exact=True), 
  hole_starts=[(0,0)])
test_tools.adaptive_demo(sol, mesh, max_iterations = 8)
