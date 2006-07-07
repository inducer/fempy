import math

import pylinear.array as num

import pymbolic

import fempy
import fempy.mesh
import fempy.geometry

import test_tools


# zero boundary
sol = pymbolic.parse("(1-x[0]**2)*(1-x[1]**2)")

# constant rhs
# sol = pymbolic.parse("x[0]**2+x[1]**2")

grid_vectors = [num.array([2,0], num.Float),
                num.array([0,2], num.Float)]
mesh = fempy.mesh.TwoDimensionalMesh(
  [fempy.mesh.ShapeSection(fempy.geometry.get_parallelogram(grid_vectors), "dirichlet")])
test_tools.adaptive_demo(sol, mesh, max_iterations = 8)
