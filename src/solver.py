import tools
import stopwatch
import element
import mesh_function
import pylinear.matrices as num
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




def getDirichletConstraints(mesh, u_d):
    constraints = {}
    for node in [node 
                 for node in mesh.dofManager()
                 if node.TrackingId == "dirichlet"]:
        constraints[node] = (u_d(node.Coordinates), [])
        print "DI!", node.Coordinates, u_d(node.Coordinates)
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
    shifted_matrix_inv_op = algo.makeCGMatrixOperator(shifted_matrix_op,
                                                      shifted_matrix_op.shape[0],
                                                      1e-15)

    op = algo.composeMatrixOperators(shifted_matrix_inv_op, m_op)

    job = stopwatch.tJob("arpack rci")
    results = algo.runArpack(op, m_op, algo.SHIFT_AND_INVERT_GENERALIZED,
                             sigma, number_of_eigenvalues, number_of_arnoldi_vectors,
                             algo.LARGEST_MAGNITUDE, tolerance, False, max_iterations)
    job.done()
  
    return zip(results.RitzValues, results.RitzVectors)




def _constraints2uc_contributions(constraints):

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

    uc_contributions = {}
    for node, (offset, lincomb_specifier) in constraints.iteritems():
        assert offset == 0
        # first, we need to make sure that we can pick out an unconstrained node.
        # so, we'll resolve all the constrained ones in lincomb_specifier down 
        # to their unconstrained definitions.
        uc_components = resolve_to_unconstrained(lincomb_specifier)

        # then, we set up the equation so that there's an unconstrained node
        # by itself
        uc_coeff, uc_node = uc_components.pop(0)
        uc_contributions[uc_node] = [(1/uc_coeff, node)] \
                                     + [(-other_coeff/uc_coeff, other_node)
                                        for other_coeff, other_node in uc_components]
    return uc_contributions

        
        

class tLaplacianEigenproblemSolver:
    def __init__(self, mesh, constrained_nodes,
                 f = None, g = lambda x: 1., typecode = num.Float):
        """Solve the eigenproblem of the Laplace operator:
        
        laplace u + f * u = g * lambda * u.
    
        For the nodes marked with the constraint id "dirichlet",
        a (boundary) condition of u = 0 is enforced.
        """

        self.Mesh = mesh
        
        dof_manager = mesh.dofManager()
        unconstrained_nodes = [node for node in dof_manager if not node in constrained_nodes]
        number_assignment = self.NumberAssignment = \
                            element.assignNodeNumbers(unconstrained_nodes)
        complete_number_assignment = self.CompleteNumberAssignment = \
                                     element.assignNodeNumbers(constrained_nodes, 
                                                               number_assignment)
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

    def massMatrix(self):
        return self.FullM

    def solve(self, sigma, constraints, number_of_eigenvalues = 20):
        dof_manager = self.Mesh.dofManager()
        full_dof_count = len(dof_manager)
        dof_count = full_dof_count - len(constraints)

        uc_contributions = _constraints2uc_contributions(constraints)

        a = num.zeros((dof_count, full_dof_count), self.FullS.typecode(), num.SparseBuildMatrix)
        for i in range(dof_count):
            a[i,i] = 1
        
        for node, lincomb_specifier in uc_contributions.iteritems():
            nonu = self.NumberAssignment[node]
            for coeff, other_node in lincomb_specifier:
                other_nonu = self.CompleteNumberAssignment[other_node]
                a[nonu, other_nonu] = coeff
        
        a_ex = num.asarray(a, a.typecode(), num.SparseExecuteMatrix)
        a_herm_ex = num.hermite(a)
        a_op = algo.makeMatrixOperator(a_ex)
        a_herm_op = algo.makeMatrixOperator(a_herm_ex)

        s_ex = num.asarray(self.FullS, self.FullS.typecode(), num.SparseExecuteMatrix)
        s_op = algo.makeMatrixOperator(s_ex)

        m_ex = num.asarray(self.FullM, self.FullM.typecode(), num.SparseExecuteMatrix)
        m_op = algo.makeMatrixOperator(m_ex)

        periodic_s_op = algo.composeMatrixOperators(
            algo.composeMatrixOperators(a_op, s_op),
            a_herm_op)

        periodic_m_op = algo.composeMatrixOperators(
            algo.composeMatrixOperators(a_op, m_op),
            a_herm_op)

        eigen_result = shiftAndInvertSymmetricEigenproblem(sigma, 
                                                           periodic_s_op, 
                                                           periodic_m_op,
                                                           number_of_eigenvalues)

        return [(eigenvalue, mesh_function.tMeshFunction(
            self.Mesh,
            self.NumberAssignment,
            mesh_function.makeCompleteVector(self.CompleteNumberAssignment,
                                             eigenvector,
                                             constraints)))
                for eigenvalue, eigenvector in eigen_result]

    def computeEigenpairResidual(self, value, mesh_func):
        s = self.FullS
        m = self.FullM

        vector = mesh_func.vector()

        if s is None or m is None:
            raise RuntimeError, "insertConstraints was never called, "+\
                  "no finalized matrices available."
        
        return tools.norm2(
            num.matrixmultiply(s, vector) 
            - value * num.matrixmultiply(m, vector))
