import math

import test_tools
import fempy
import fempy.mesh
import fempy.geometry
import fempy.expression as expression

outer_radius = 5
inner_radius = 0.1
outer_value = 0
inner_value = 10

c1 = (outer_value - inner_value)/(math.log(outer_radius/inner_radius))
c2 = inner_value - c1 * math.log(inner_radius)

expr_r_squared = ("+", ("**", ("variable", "0"), 2), ("**", ("variable", "1"), 2))
expr_r = ("**", expr_r_squared, 0.5)
sol = ("+", ("*", c1, ("log", expr_r)), c2)

mesh = fempy.mesh.tTwoDimensionalMesh(
  fempy.geometry.getAnnulusGeometry(5, 0.1, use_exact = True), 
  hole_starts = [(0,0)])
test_tools.adaptiveDemo(sol, mesh, max_iterations = 6)
