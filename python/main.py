#! /usr/bin/python

import Numeric as num
import spmatrix as sparse




class tMatrixBuilder:
  def __init__( self, size_hint = 1000 ):
    self.Matrix = None
    self.Unknonws = { }

  def obtainUnknown( self, identifier ):
    if identifier in self.Unknowns:
      return self.Unknowns[ identifier ]
    else:
      new_unknown_id = len( self.Unknowns )
      self.Unknowns[ identifier ] = new_unknown_id
      return new_unknown_id
    
  def beginBuildingMatrix( self ):
    pass
    self.Matrix = sparse.ll_mat( len( self.Unknowns ), self.UnknownCount )

  def addToMatrix( self, node_ids, little_matrix ):
    pass
    self.Matrix.update_add_mask( little_matrix, node_ids, node_ids )




class tSymmetricSparseMatrixBuilder( tMatrixBuilder ):
  pass




class tGeneralSparseMatrixBuilder( tMatrixBuilder ):
  def beginBuildingMatrix( self ):
    self.Matrix = sparse.ll_mat( len( self.Unknowns ), self.UnknownCount )

  def addToMatrix( self, node_ids, little_matrix ):
    self.Matrix.update_add_mask( little_matrix, node_ids, node_ids )




def tVectorBuilder( tMatrixBuilder ):
  def beginBuildingMatrix( self ):
    self.Matrix = num.zeros( (self.UnknownCount,) )

  def addToMatrix( self, node_ids, little_matrix ):
    for i in range( 0, len( node_ids ) ):
      little_matrix[ node_ids[ i ] ] += little_matrix[ i ]




# geometry classes ------------------------------------------------------------
class tNode:
  def __init__( self, coordinates ):
    self.Coordinates = coordinates

  def coordinates( self ):
    return self.Coordinates




# finite element classes ------------------------------------------------------
class tFiniteElementError(Exception):
  def __init__(self, value):
    self.value = value
  def __str__(self):
    return repr(self.value)




class tFiniteElement:
  def __init__( self, nodes ):
    self.Nodes = nodes

  def allocateUnknowns( self, builder ):
    pass

  def addVolumeIntegralOverUxSquared( self, builder ):
    pass

  def addVolumeIntegralOverUySquared( self, builder ):
    pass
  
  def addVolumeIntegralOverUSquared( self, builder ):
    pass

  def addVolumeIntegralOverU( self, builder ):
    pass





class tTriangularFiniteElement( tFiniteElement ):
  def __init__( self, nodes ):
    tFiniteElement.__init__( self, nodes )
    self.X = map( self.Nodes, lambda node: node.coordinates()[ 0 ] )
    self.Y = map( self.Nodes, lambda node: node.coordinates()[ 1 ] )

    self.Area = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0])

  def area( self ):
    return self.Area

  def barycenter( self ):
    return 





class tTwoDimensionalLinearTriangularFiniteElement( tFiniteElement ):
  # cf. \cite{50linesofmatlab} 
  # also see related axiom calculation 50_lines_of_matlab_matrix_m.tm

  def __init__( self, nodes ):
    tTriangularFiniteElement.__init__( self, nodes )
    x = self.X
    y = self.Y

    self.InvCoordinateMatrix = la.inverse( num.concatenate(
       num.ones( (1,3) ), num.array( [x,y] ) ) )

  def allocateUnknowns( self, builder ):
    node_ids = map( builder.obtainUnknown, self.Nodes )

  def addVolumeIntegralOverUxSquared( self, builder ):
    node_ids = map( builder.obtainUnknown, self.Nodes )
    gx = num.matrixmultiply( self.InvCoordinateMatrix, num.array(
      [ [ 0 0 ], [ 0 0 ], [ 0 1] ] ) )
    influence_matrix = self.Area * 0.5 * g * transpose( g )
    builder.addToMatrix( node_ids, influence_matrix )

  def addVolumeIntegralOverUySquared( self, builder ):
    node_ids = map( builder.obtainUnknown, self.Nodes )
    gx = num.matrixmultiply( self.InvCoordinateMatrix, num.array(
      [ [ 0 0 ], [ 1 0 ], [ 0 0] ] ) )
    influence_matrix = self.Area * 0.5 * g * transpose( g )
    builder.addToMatrix( node_ids, influence_matrix )
  
  def addVolumeIntegralOverUSquared( self, builder ):
    raise tFiniteElementError, "NYI: can't add u^2 yet"

  def addVolumeIntegralOverU( self, builder ):
    pass




# equation solvers ------------------------------------------------------------





# visualization ---------------------------------------------------------------
# geometry builder ------------------------------------------------------------





# driver ----------------------------------------------------------------------
