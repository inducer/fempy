import math
import cmath

# Numerics imports ------------------------------------------------------------
import pylinear.matrices as num
import pylinear.linear_algebra as la

# fempy -----------------------------------------------------------------------
import fempy.mesh
import fempy.geometry
import fempy.stopwatch
import fempy.solver
import fempy.tools as tools

import test_tools




def unitCellDemo(mesh, epsilon, sigma, k):
  job = fempy.stopwatch.tJob("geometry")
  mesh.generate()
  job.done()

  boundaries = [[-0.5, 0.5], [-0.5,0.5]]
  def isOnSameBoundary(point1, point2):
    n, = point1.shape
    for i in range(n):
      if point1[i] in boundaries[i] and point2[i] == point1[i]:
        return True
    return False

  job = fempy.stopwatch.tJob("periodicity")
  thresh = math.pi / 8.
  bnodes = mesh.boundaryNodes()

  pairs = []
  for index, node in zip(range(len(bnodes)), bnodes):
    without_this = bnodes[index+1:]

    nodes_with_distances = tools.decorate(
      lambda other_node: tools.angleBetweenVectors(node.coordinates()-other_node.coordinates(), k),
      without_this)
    nodes_with_distances.sort(lambda (node1, dist1), (node2, dist2): cmp(dist1, dist2))
    
    taken = 0 
    for other_node, dist in nodes_with_distances:
      if dist > thresh and taken >= 5:
        break
      pairs.append((node, other_node, dist))

  pairs.sort(lambda (node1a, node1b, dist1), (node2a, node2b, dist2): cmp(dist1, dist2))

  unconnected = tools.tReference(len(bnodes))
  periodic_nodes = []
  connections = {}
  
  def connect(node_a,node_b):
    if node_a not in connections:
      unconnected.set(unconnected.get() - 1)
      connections[node_a] = [node_b]
    else:
      connections[node_a].append(node_b)
    if node_b not in connections:
      unconnected.set(unconnected.get() - 1)
      connections[node_b] = [node_a]
    else:
      connections[node_b].append(node_a)

  pairs_index = 0
  while unconnected.get() > 0 and pairs_index < len(pairs):
    node_a, node_b, dist = pairs[pairs_index]
    if (not node_a in connections or \
       not node_b in connections) and \
       dist < thresh:
      connect(node_a, node_b)

      k_dist = num.innerproduct(node_a.coordinates()-node_b.coordinates(),k)
      periodic_nodes.append((node_a, node_b, cmath.exp(1j * k_dist)))
    pairs_index += 1

  connections_file = file(",,connections.data", "w")
  for a,b,factor in periodic_nodes:
    connections_file.write("%f\t%f\n" % (a.coordinates()[0], a.coordinates()[1]))
    connections_file.write("%f\t%f\n" % (b.coordinates()[0], b.coordinates()[1]))
    connections_file.write("\n")
  job.done()

  results = fempy.solver.solveLaplaceEigenproblem(sigma, mesh, 
                                                  [], periodic_nodes, 
                                                  g = epsilon,
                                                  typecode = num.Complex)

  for val, vec in zip(results.RitzValues, results.RitzVectors):
    print val
    test_tools.visualize(mesh, num.real(vec))
    raw_input()




def runEigenDemo():
  def needsRefinement( vert_origin, vert_destination, vert_apex, area ):
    bary_x = ( vert_origin.x() + vert_destination.x() + vert_apex.x() ) / 3
    bary_y = ( vert_origin.y() + vert_destination.y() + vert_apex.y() ) / 3
    
    dist_center = math.sqrt( bary_x**2 + bary_y**2 )
    if dist_center < 0.4:
      return area >= 1e-3
    else:
     return False

  mesh = fempy.mesh.tTwoDimensionalMesh(
    fempy.geometry.getUnitCellGeometry(edge_length = 1, 
                                 inner_factor = 0.3,
                                 use_exact = False),
    refinement_func = needsRefinement)

  def epsilon(x):
    if tools.norm2(x) < 0.3:
      return 11
    else:
      return 0

  sigma = 0.9
  
  unitCellDemo(mesh, epsilon, sigma, num.array([1,0.5]))




runEigenDemo()
