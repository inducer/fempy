import math
import tools
import stopwatch
import element
import mesh_function
import visualization
import pylinear.matrices as num
import pylinear.matrix_tools as mtools
import pylinear.algorithms as algo




def solveSPDSystem(matrix_op, rhs, start_vector = None):
    h,w = matrix_op.shape
    if start_vector is None:
        x = num.zeros((w,), num.Float)
    else:
        x = start_vector[:]

    job = stopwatch.tJob("solve")
    matrix_inv_op = algo.makeCGMatrixOperator(matrix_op, h * 2)
    matrix_inv_op.apply(rhs, x)
    job.done()
    print "  iterations:", matrix_inv_op.last_iteration_count

    residual = x.copy()
    matrix_op.apply(x, residual)
    residual -= rhs

    print "  absolute residual: ", tools.norm2(residual)
    return x




def getNodesWithTrackingId(mesh, id):
    return [node 
            for node in mesh.dofManager()
            if node.TrackingId == id]




def getDirichletConstraints(mesh, u_d):
    constraints = {}
    for node in [node 
                 for node in mesh.dofManager()
                 if node.TrackingId == "dirichlet"]:
        constraints[node] = (u_d(node.Coordinates), [])

    return constraints




def buildStiffnessMatrix(mesh, number_assignment, typecode = num.Float):
    dof_count = len(number_assignment)
    s = num.zeros((dof_count, dof_count), typecode, num.SparseBuildMatrix)

    job = stopwatch.tJob("stiffness matrix")
    for el in mesh.elements():
        nonus = [number_assignment[node] for node in el.nodes()]
        s.addScatteredSymmetric(
            nonus, 
            num.asarray(el.getVolumeIntegralsOverDifferentiatedFormFunctions(),
                        typecode))

    job.done()

    return s




def buildMassMatrix(mesh, number_assignment, weight_function, typecode = num.Float):
    dof_count = len(number_assignment)
    m = num.zeros((dof_count, dof_count), typecode, num.SparseBuildMatrix)

    job = stopwatch.tJob("mass matrix")
    for el in mesh.elements():
        nonus = [number_assignment[node] for node in el.nodes()]
        m.addScatteredSymmetric(
            nonus, 
            num.asarray(el.getVolumeIntegralsOverFormFunctions(weight_function),
                        typecode))

    job.done()

    return m




def solvePoisson(mesh, f, node_constraints, start_vector = None):
    """Solve the Poisson equation

    laplace u = f
    with u = u_d in the dirichlet_nodes.

    node_constraints is a dictionary that associates nodes with 
    constraints. A constraint is a datastructure of the following shape:
    
    intercept, [(coefficient_1, node_2), (coefficient_2, index_2), ...]

    which forces the node to be 

    intercept + coefficient_1 * getValue(node_1) + ...
    """

    dof_manager = mesh.dofManager()
    unconstrained_nodes = [node for node in dof_manager if not node in node_constraints]
    constrained_nodes = node_constraints.keys()
    number_assignment = element.assignNodeNumbers(unconstrained_nodes)
    complete_number_assignment = element.assignNodeNumbers(constrained_nodes, 
                                                           number_assignment)
    dof_count = len(unconstrained_nodes)

    full_s = buildStiffnessMatrix(mesh, complete_number_assignment)
    s = full_s[:dof_count,:dof_count]

    b = num.zeros((dof_count,), num.Float)

    job = stopwatch.tJob("rhs")
    for el in mesh.elements():
        this_el_b = el.getVolumeIntegralsOverFormFunction(
            lambda x,formfunc_value: f(x) * formfunc_value)
        for node, v in zip(el.nodes(), this_el_b):
            try:
                b[number_assignment[node]] += v
            except KeyError:
                pass
    job.done()

    job = stopwatch.tJob("bcs: dirichlet")
    for node, (boundary_value, other_nodes) in node_constraints.iteritems():
        assert len(other_nodes) == 0
        if boundary_value != 0:
            nonu = complete_number_assignment[node]
            b += boundary_value * full_s[:dof_count, nonu]
    job.done()

    #visualization.writeGnuplotSparsityPattern(",,s.data", s)

    compiled_s = num.asarray(s, s.typecode(), num.SparseExecuteMatrix)
    s_op = algo.makeMatrixOperator(compiled_s)
    
    complete_vec = mesh_function.makeCompleteVector(
        complete_number_assignment,
        solveSPDSystem(s_op, -b),
        node_constraints)
    return mesh_function.tMeshFunction(
        mesh, complete_number_assignment, complete_vec)




def shiftAndInvertEigenproblem(sigma, s, m, 
                               number_of_eigenvalues = 20, 
                               tolerance = 1e-10,
                               max_iterations = 0):

    number_of_arnoldi_vectors = 2 * number_of_eigenvalues
    compiled_m = num.asarray(m, m.typecode(), num.SparseExecuteMatrix)

    m_op = algo.makeMatrixOperator(compiled_m)

    job = stopwatch.tJob("shift matrix")
    shifted_matrix = num.asarray(s - sigma * m, s.typecode(), num.SparseExecuteMatrix)
    shifted_matrix_invop = algo.makeUMFPACKMatrixOperator(shifted_matrix)
    job.done()

    op = algo.composeMatrixOperators(shifted_matrix_invop, m_op)

    job = stopwatch.tJob("arpack rci")
    results = algo.runArpack(op, m_op, algo.SHIFT_AND_INVERT_GENERALIZED,
                             sigma, number_of_eigenvalues, number_of_arnoldi_vectors,
                             algo.LARGEST_MAGNITUDE, tolerance, False, max_iterations)
    job.done()
  
    return zip(results.RitzValues, results.RitzVectors)




def shiftAndInvertSymmetricEigenproblem(sigma, s_op, m_op, 
                               number_of_eigenvalues = 20, 
                               tolerance = 1e-10,
                               max_iterations = 0):

    number_of_arnoldi_vectors = 2 * number_of_eigenvalues

    neg_sigma_m_op = algo.composeMatrixOperators(
        algo.makeScalarMultiplicationMatrixOperator(-sigma, 
                                                    s_op.shape[0], 
                                                    s_op.typecode()),
        m_op)
    shifted_matrix_op = algo.addMatrixOperators(s_op, neg_sigma_m_op)
    shifted_matrix_inv_op = algo.makeBiCGSTABMatrixOperator(
        shifted_matrix_op, shifted_matrix_op.shape[0]*10, tolerance)
    #shifted_matrix_inv_op.debug_level = 2

    op = algo.composeMatrixOperators(shifted_matrix_inv_op, m_op)

    job = stopwatch.tJob("arpack rci")
    results = algo.runArpack(op, m_op, algo.SHIFT_AND_INVERT_GENERALIZED,
                             sigma, number_of_eigenvalues, number_of_arnoldi_vectors,
                             algo.LARGEST_MAGNITUDE, tolerance, False, max_iterations)
    job.done()
  
    return zip(results.RitzValues, results.RitzVectors)




def resolveConstraints(constraints):
    """Resolves the right hand sides of the constraint equations 
    to only unconstrained nodes.
    """

    def resolve_to_unconstrained(lincomb_specifier):
        uc_components = [(other_coeff, other_node)
                         for other_coeff, other_node in lincomb_specifier
                         if other_node not in constraints]
        constr_components = [(other_coeff, other_node)
                         for other_coeff, other_node in lincomb_specifier
                         if other_node in constraints]

        while constr_components:
            other_coeff, other_node = constr_components.pop(0)
            other_offset, other_lincomb_specifier = constraints[other_node]
            for other2_coeff, other2_node in other_lincomb_specifier:
                if other2_node in constraints:
                    constr_components.append((other_coeff * other2_coeff, other2_node))
                else:
                    uc_components.append((other_coeff * other2_coeff, other2_node))
        return uc_components

    new_constraints = {}
    for node, (offset, lincomb_specifier) in constraints.iteritems():
        assert offset == 0
        #print "resolving", node.Coordinates, "out of",
        #for coeff, other_node in lincomb_specifier:
            #print other_node.Coordinates,
        #print

        new_constraints[node] = resolve_to_unconstrained(lincomb_specifier)

    return new_constraints




def rayleighQuotient(s, m, vec):
    sp = mtools.sp
    mm = num.matrixmultiply
    return sp(vec, mm(s, vec))/sp(vec, mm(m, vec))




class tLaplacianEigenproblemSolver:
    def __init__(self, mesh, constrained_nodes,
                 f = None, g = lambda x: 1., typecode = num.Float,
                 given_number_assignment = None):
        """Solve the eigenproblem of the Laplace operator:
        
        laplace u + f * u = g * lambda * u.
        """

        self.Mesh = mesh
        
        dof_manager = mesh.dofManager()
        unconstrained_nodes = [node for node in dof_manager if node not in constrained_nodes]
        if given_number_assignment is None:
            # construct new node number assignment
            number_assignment = self.NumberAssignment = \
                                element.assignNodeNumbers(unconstrained_nodes)
            complete_number_assignment = self.CompleteNumberAssignment = \
                                         element.assignNodeNumbers(constrained_nodes, 
                                                                   number_assignment)
        else:
            # make sure that the given node number assignment conforms to our convention
            assert max(given_number_assignment.values()) == len(mesh.dofManager()) - 1

            complete_number_assignment = self.CompleteNumberAssignment = \
                                         given_number_assignment
            number_assignment = self.NumberAssignment = {}
            for node in unconstrained_nodes:
                assert given_number_assignment[node] < len(unconstrained_nodes)
                number_assignment[node] = given_number_assignment[node]
                                
        dof_count = len(unconstrained_nodes)

        self.FullS = buildStiffnessMatrix(mesh, complete_number_assignment, typecode)
        self.FullM = buildMassMatrix(mesh, complete_number_assignment, g, typecode)

        if f is not None:
            job = stopwatch.tJob("f")
    
            for el in mesh.elements():
                nonus = [number_assignment[node] for node in el.nodes()]
                self.FullS.addScatteredSymmetric(
                    nonus, 
                    num.asarray(el.getVolumeIntegralsOverFormFunctions(f), typecode))
            job.done()

        self.CurrentConstraints = None
        self.ConstrainedSOp = self.ConstrainedMOp = None

    def stiffnessMatrix(self):
        return num.asarray(self.FullS, self.FullS.typecode(), num.SparseExecuteMatrix)

    def massMatrix(self):
        return num.asarray(self.FullM, self.FullM.typecode(), num.SparseExecuteMatrix)

    def nodeNumberAssignment(self):
        return self.CompleteNumberAssignment

    def currentConstraints(self):
        return self.CurrentConstraints

    def setupConstraints(self, constraints):
        dof_manager = self.Mesh.dofManager()
        full_dof_count = len(dof_manager)
        dof_count = full_dof_count - len(constraints)

        resolved_constraints = resolveConstraints(constraints)

        a = num.zeros((dof_count, full_dof_count), self.FullS.typecode(), 
                      num.SparseBuildMatrix)
        for i in range(dof_count):
            a[i,i] = 1

        for c_node, lincomb_specifier in resolved_constraints.iteritems():
            c_nonu = self.CompleteNumberAssignment[c_node]
            assert c_node not in self.NumberAssignment
            for coeff, uc_node in lincomb_specifier:
                uc_nonu = self.NumberAssignment[uc_node]
                a[uc_nonu, c_nonu] += coeff
        
        a_ex = num.asarray(a, a.typecode(), num.SparseExecuteMatrix)
        a_herm_ex = num.hermite(a_ex)
        a_op = algo.makeMatrixOperator(a_ex)
        a_herm_op = algo.makeMatrixOperator(a_herm_ex)

        #mm = num.matrixmultiply
        #print "A"
        #mtools.printComplexMatrixInGrid(a)
        #print "S"
        #mtools.printComplexMatrixInGrid(self.FullS)
        #print "M"
        #mtools.printComplexMatrixInGrid(self.FullM)
        #print "ASA_H"
        #mtools.printComplexMatrixInGrid(mm(a, mm(self.FullS, num.hermite(a))))
        #print "AMA_H"
        #mtools.printComplexMatrixInGrid(mm(a, mm(self.FullM, num.hermite(a))))

        s_ex = num.asarray(self.FullS, self.FullS.typecode(), num.SparseExecuteMatrix)
        s_op = algo.makeMatrixOperator(s_ex)

        m_ex = num.asarray(self.FullM, self.FullM.typecode(), num.SparseExecuteMatrix)
        m_op = algo.makeMatrixOperator(m_ex)

        self.CurrentConstraints = constraints
        self.ConstrainedSOp = algo.composeMatrixOperators(
            algo.composeMatrixOperators(a_op, s_op),
            a_herm_op)
        self.ConstrainedMOp = algo.composeMatrixOperators(
            algo.composeMatrixOperators(a_op, m_op),
            a_herm_op)
        return a

    def solve(self, sigma, number_of_eigenvalues = 20, tolerance = 1e-10,
              warning_threshold = 10):
        if self.CurrentConstraints is None:
            raise RuntimeError, "need to set up constraints before solving"

        eigen_result = shiftAndInvertSymmetricEigenproblem(sigma, 
                                                           self.ConstrainedSOp, 
                                                           self.ConstrainedMOp,
                                                           number_of_eigenvalues,
                                                           tolerance = tolerance)

        result = [(eigenvalue, mesh_function.tMeshFunction(
            self.Mesh,
            self.CompleteNumberAssignment,
            mesh_function.makeCompleteVector(self.CompleteNumberAssignment,
                                             num.conjugate(eigenvector),
                                             self.CurrentConstraints)))
                for eigenvalue, eigenvector in eigen_result]

        # perform some invariant checking
        if warning_threshold is not None:
            for evalue, mf in result:
                if self.computeEigenpairResidual(evalue, mf) > warning_threshold * tolerance:
                    print "PRECSION WARNING for ARPACK output."

        return result


    def computeEigenpairResidual(self, eigenvalue, mesh_func):
        if self.CurrentConstraints is None:
            raise RuntimeError, "need to set up constraints before verifying"

        dof_count = len(self.Mesh.dofManager())- len(self.CurrentConstraints)
        vector = num.conjugate(mesh_func.vector()[:dof_count])

        def m_sp(x,y): return mtools.sp(x, algo.applyMatrixOperator(
            self.ConstrainedMOp, y))
        def m_norm(x): 
            return math.sqrt(abs(m_sp(x,x)))

        return m_norm(
            algo.applyMatrixOperator(self.ConstrainedSOp, vector)
            - eigenvalue * algo.applyMatrixOperator(self.ConstrainedMOp, vector))
