import pylinear.matrices as num




class tVectorBuilder:
  def matrix(self):
    """The result of this function is, in general, of unspecified type.
    However, the returned object will at least support the following
    subset of the numarray interface:

      m.shape
      subscripting and slicing
    """
    pass

  def addScattered(self, small_matrix, small_matrix_rows):
    pass




class tMatrixBuilder:
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
    pass

  def addScattered(self, small_matrix, small_matrix_rows, small_matrix_columns):
    pass

  def addScatteredSymmetric(self, small_matrix, small_matrix_rows):
    self.addScattered(small_matrix, small_matrix_rows, small_matrix_rows)





class tPyLinearMatrixBuilder(tMatrixBuilder):
  def __init__(self):
    self.Matrix = None

  def matrix(self):
    return self.Matrix

  def forceIdentityMap(self, dof_number):
    self.Matrix[dof_number] = 0
    self.Matrix[:,dof_number] = 0
    self.Matrix[dof_number, dof_number] = 1.

  def column(self, i):
    return self.Matrix[:,i]

  def addScattered(self, small_matrix, small_matrix_rows, small_matrix_columns):
    self.Matrix.addScattered(small_matrix_rows, small_matrix_columns, 
                             num.asarray(small_matrix, self.Matrix.typecode()))

  def addScatteredSymmetric(self, small_matrix, small_matrix_rows):
    self.Matrix.addScatteredSymmetric(small_matrix_rows, 
                                      num.asarray(small_matrix, self.Matrix.typecode()))




class tSymmetricSparseMatrixBuilder(tPyLinearMatrixBuilder):
  def __init__(self, size, typecode):
    tPyLinearMatrixBuilder.__init__(self)
    self.Matrix = num.zeros((size, size), typecode, num.SparseBuildMatrix)

  def addScattered(self, small_matrix, small_matrix_rows, small_matrix_columns):
    raise RuntimeError, "cannot addScattered to a symmetric matrix"




class tGeneralSparseMatrixBuilder(tPyLinearMatrixBuilder):
  def __init__(self, shape, typecode):
    tPyLinearMatrixBuilder.__init__(self)
    self.Matrix = num.zeros(shape, typecode, num.SparseBuildMatrix)




class tDenseMatrixBuilder(tPyLinearMatrixBuilder):
  def __init__(self, shape, typecode):
    tPyLinearMatrixBuilder.__init__(self)
    self.Matrix = num.zeros(shape, typecode)




class tDenseVectorBuilder(tVectorBuilder):
  def __init__(self, size, typecode):
    self.Matrix = num.zeros((size,), typecode)

  def matrix(self):
    return self.Matrix

  def addScattered(self, small_matrix, small_matrix_rows):
    for i in range(0, len(small_matrix_rows)):
      if small_matrix[ i ] != 0:
	self.Matrix[ small_matrix_rows[ i ] ] += small_matrix[ i ]

