import Numeric as num
import LinearAlgebra as la
import spmatrix as sparse




class tDOFManager:
  def __init__(self):
    self.IdentifierToNumber = { }
    self.NumberToIdentifier = [ ]

  def getDegreeOfFreedomNumber(self, identifier):
    if identifier in self.IdentifierToNumber:
      return self.IdentifierToNumber[ identifier ]
    else:
      new_dof_id = len(self.NumberToIdentifier)
      self.IdentifierToNumber[ identifier ] = new_dof_id
      self.NumberToIdentifier.append(identifier)
      return new_dof_id

  def getDegreeOfFreedomIdentifier(self, number):
    return self.NumberToIdentifier[ number ]

  def countDegreesOfFreedom(self):
    return len(self.NumberToIdentifier)




class tMatrixBuilder:
  def add(self, small_matrix, small_matrix_rows, small_matrix_columns = None):
    pass

  def matrix(self):
    """The result of this function is, in general, of unspecified type.
    However, the returned object will at least support the following
    subset of the numarray interface:

      m.shape
      subscripting and slicing
    """
    pass

  def column(self, i):
    """Returns the i'th column of the built matrix as a copied dense array."""
    pass

  def forceIdentityMap(self, dof_number):
    mat = self.Matrix
    h,w = mat.shape

    for i in range(0, w):
      mat[ dof_number, i ] = 0
    for i in range(0, h):
      mat[ i, dof_number ] = 0
    mat[ dof_number, dof_number ] = 1.




class tSymmetricSparseMatrixBuilder(tMatrixBuilder):
  def __init__(self, size):
    self.Matrix = sparse.ll_mat_sym(size)

  def matrix(self):
    return self.Matrix

  def forceIdentityMap(self, dof_number):
    mat = self.Matrix
    h,w = mat.shape

    # FIXME: optimize with slicing syntax
    for i in range(0, dof_number):
      mat[ dof_number, i ] = 0
    for i in range(dof_number + 1, h):
      mat[ i, dof_number ] = 0
    mat[ dof_number, dof_number ] = 1.

  def column(self, i):
    h,w = self.Matrix.shape
    col = num.zeros((h,), num.Float)
    for j in range(0, i):
      col[ j ] = self.Matrix[ i,j ]
    for j in range(i, h):
      col[ j ] = self.Matrix[ j,i ]
    return col

  def add(self, small_matrix, small_matrix_rows, small_matrix_columns = None):
    self.Matrix.update_add_mask_sym(
	small_matrix, 
	num.array(small_matrix_rows),
	num.ones((len(small_matrix_rows),)))




class tGeneralSparseMatrixBuilder(tMatrixBuilder):
  def __init__(self, height, width):
    self.Matrix = sparse.ll_mat(height, width)

  def matrix(self):
    return self.Matrix

  def add(self, small_matrix, small_matrix_rows, small_matrix_columns = None):
    if small_matrix_columns is None:
      small_matrix_columns = small_matrix_rows
    self.Matrix.update_add_mask(small_matrix, num.array(small_matrix_rows), num.array(small_matrix_columns))




class tDenseVectorBuilder(tMatrixBuilder):
  def __init__(self, size):
    self.Matrix = num.zeros((size,), num.Float)

  def matrix(self):
    return self.Matrix

  def add(self, small_matrix, small_matrix_rows, small_matrix_columns = None):
    for i in range(0, len(small_matrix_rows)):
      if small_matrix[ i ] != 0:
	self.Matrix[ small_matrix_rows[ i ] ] += small_matrix[ i ]




class tDenseMatrixBuilder(tMatrixBuilder):
  def __init__(self, height, width):
    self.Matrix = num.zeros((height,width), num.Float)

  def matrix(self):
    return self.Matrix

  def add(self, small_matrix, small_matrix_rows, small_matrix_columns = None):
    h,w = self.Matrix.shape
    if small_matrix_columns is None:
      small_matrix_columns = small_matrix_rows

    for i in range(0, len(small_matrix_rows)):
      for j in range(0, len(small_matrix_columns)):
	self.Matrix[ small_matrix_rows[ i ], small_matrix_columns[ j ] ] += small_matrix[ i ][ j ]

  def forceIdentityMap(self, dof_number):
    mat = self.Matrix
    mat[ :,dof_number ] = 0
    mat[ dof_number ] = 0
    mat[ dof_number, dof_number ] = 1.





