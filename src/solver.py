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
    return constraints




def _buildStiffnessMatrix(mesh, dof_count, number_assignment, typecode = num.Float):
    s = num.zeros((dof_count, dof_count), typecode, num.SparseBuildMatrix)

    contrib_of_constrained_nodes = tools.tDictionaryWithDefault(
        lambda x: tools.tSparseVector())
    
    job = stopwatch.tJob("stiffness matrix")
    for el in mesh.elements():
        nonus = [number_assignment[node] for node in el.nodes()]
        contribs = num.asarray(el.getVolumeIntegralsOverDifferentiatedFormFunctions(),
                               typecode)
        s.addScatteredSymmetric(nonus, contribs)

        constrained_indices = [i for i in range(len(nonus)) if nonus[i] < 0]
        if constrained_indices:
            nodes = el.nodes()
            unconstrained_indices = [i for i in range(len(nonus)) if nonus[i] >= 0]
            for i in constrained_indices:
                for j in unconstrained_indices:
                    contrib_of_constrained_nodes[nodes[i]][nonus[j]] += contribs[i,j]

    job.done()

    return s, contrib_of_constrained_nodes




def _buildMassMatrix(mesh, dof_count, number_assignment, weight_function, typecode = num.Float):
    m = num.zeros((dof_count, dof_count), typecode, num.SparseBuildMatrix)

    contrib_of_constrained_nodes = tools.tDictionaryWithDefault(
        lambda x: tools.tSparseVector())
    
    job = stopwatch.tJob("mass matrix")
    for el in mesh.elements():
        nonus = [number_assignment[node] for node in el.nodes()]
        contribs = num.asarray(el.getVolumeIntegralsOverFormFunctions(weight_function),
                               typecode)
        m.addScatteredSymmetric(nonus, contribs)

        constrained_indices = [i for i in range(len(nonus)) if nonus[i] < 0]
        if constrained_indices:
            nodes = el.nodes()
            unconstrained_indices = [i for i in range(len(nonus)) if nonus[i] >= 0]
            for i in constrained_indices:
                for j in unconstrained_indices:
                    contrib_of_constrained_nodes[nodes[i]][nonus[j]] += contribs[i,j]

    job.done()

    return m, contrib_of_constrained_nodes




def solvePoisson(mesh, f, node_constraints, start_vector = None):
    """Solve the Poisson equation

    laplace u = f
    with u = u_d in the dirichlet_nodes.

    node_constraints is a dictionary that associates nodes with 
    constraints. A constraint is a datastructure of the following shape:
    
    intercept, [(coefficient_1, node_2), (coefficient_2, index_2), ...]

    which is evaluated to

    intercept + coefficient_1 * getValue(node_1) + ...
    """

    dof_manager = mesh.dofManager()
    unconstrained_nodes = [node for node in dof_manager if not node in node_constraints]
    constrained_nodes = node_constraints.keys()
    number_assignment = element.assignNodeNumbers(unconstrained_nodes)
    complete_number_assignment = element.assignNodeNumbers(constrained_nodes, 
                                                           number_assignment, -1, -1)
    dof_count = len(unconstrained_nodes)

    s, contrib_of_constrained_nodes = _buildStiffnessMatrix(mesh, 
                                                            dof_count,
                                                            complete_number_assignment)

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
            contrib_of_constrained_nodes[node].addTo(b, boundary_value)
    job.done()

    #visualization.writeGnuplotSparsityPattern(",,s.data", s)

    compiled_s = num.asarray(s, s.typecode(), num.SparseExecuteMatrix)
    s_op = algo.makeMatrixOperator(compiled_s)
    
    return mesh_function.tMeshFunction(mesh, number_assignment, 
                                       solveSPDSystem(s_op, -b),
                                       node_constraints)




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




class tLaplacianEigenproblemSolver:
    def __init__(self, mesh, f = None, g = lambda x: 1., typecode = num.Float):
        """Solve the eigenproblem of the Laplace operator:
        
        laplace u + f * u = g * lambda * u.
    
        For the nodes marked with the constraint id "dirichlet",
        a (boundary) condition of u = 0 is enforced.
        """

        self.Mesh = mesh
        
        dof_manager = mesh.dofManager()
        number_assignment = self.NumberAssignment = \
                            element.assignNodeNumbers(dof_manager.unconstrainedNodes())
        complete_number_assignment = element.assignNodeNumbers(dof_manager.constrainedNodes(), 
                                                               number_assignment, 
                                                               -1, -1)
        dof_count = len(dof_manager.unconstrainedNodes())

        self.S, self.SContributionsFromConstrainedNodes = _buildStiffnessMatrix(mesh, complete_number_assignment, typecode)
        self.M, self.MContributionsFromConstrainedNodes = _buildMassMatrix(mesh, complete_number_assignment, g, typecode)
        self.SWithConstraints = None
        self.MWithConstraints = None

        if f is not None:
            job = stopwatch.tJob("rhs")
    
            for el in mesh.elements():
                nonus = [number_assignment[node] for node in el.nodes()]
                contribs = num.asarray(el.getVolumeIntegralsOverFormFunctions(f),
                                       typecode)
                s.addScatteredSymmetric(nonus, contribs)
                
                constrained_indices = [i for i in range(len(nonus)) if nonus[i] < 0]
                if constrained_indices:
                    nodes = el.nodes()
                    unconstrained_indices = [i for i in range(len(nonus)) if nonus[i] >= 0]
                    for i in constrained_indices:
                        for j in unconstrained_indices:
                            self.SContributionsFromConstrainedNodes[nodes[i]][nonus[j]] += contribs[i,j]
            job.done()

    def insertConstraints(self):
        """Update the matrices to reflect the constraints of the mesh nodes.
        May be called multiple times without harm.
        """

        s = self.SWithConstraints = self.S.copy()
        m = self.MWithConstraints = self.M.copy()

        dof_manager = self.Mesh.dofManager()

        s_contrib = self.SContributionsFromConstrainedNodes
        m_contrib = self.MContributionsFromConstrainedNodes

        nuass = self.NumberAssignment

        def addContributions(independent_node, coeff, dependent_node):
            if independent_node.Constraint is None:
                indep_node_number = self.NumberAssignment[independent_node]
                s_contrib[dependent_node].addToMatrixColumn(s, indep_node_number,
                                                             coeff)
                s_contrib[dependent_node].addToMatrixRow(s, indep_node_number,
                                                          1/coeff)
                m_contrib[dependent_node].addToMatrixColumn(m, indep_node_number,
                                                             coeff)
                m_contrib[dependent_node].addToMatrixRow(m, indep_node_number,
                                                          1/coeff)
            else:
                # Ok, we've seen that our "independent node" is really
                # just another dependent node. We'll call it an "intermediate_node".
                # Note that the parameter "independent_node" is overwritten
                # subsequently.
                intermediate_node = independent_node

                # assert no constant offset
                assert intermediate_node.Constraint[0] == 0

                for sub_coeff, independent_node in intermediate_node.Constraint[1]:
                    addContributions(independent_node, coeff * sub_coeff, dependent_node)

        job = stopwatch.tJob("linear combination constraints")
        for dependent_node in [node for node 
                               in self.Mesh.dofManager() 
                               if node.Constraint is not None]:
            # assert no constant offset
            assert dependent_node.Constraint[0] == 0

            for coeff, independent_node in dependent_node.Constraint[1]:
                addContributions(independent_node, coeff, dependent_node)
        job.done()

    def massMatrix(self):
        if self.MWithConstraints is None:
            raise RuntimeError, "updateConstraints was never called, "+\
                  "no finalized matrices available."
        else:
            return self.MWithConstraints

    def solve(self, sigma, number_of_eigenvalues = 20):
        s = self.SWithConstraints
        m = self.MWithConstraints

        if s is None or m is None:
            raise RuntimeError, "insertConstraints was never called, "+\
                  "no finalized matrices available."
        
        eigen_result = shiftAndInvertEigenproblem(sigma, s, m,
                                                  number_of_eigenvalues)

        return [(eigenvalue, mesh_function.tMeshFunction(self.Mesh,
                                                        self.NumberAssignment,
                                                        eigenmode))
                for eigenvalue, eigenmode in eigen_result]

    def computeEigenpairResidual(self, value, vector):
        s = self.SWithConstraints
        m = self.MWithConstraints

        if s is None or m is None:
            raise RuntimeError, "insertConstraints was never called, "+\
                  "no finalized matrices available."
        
        return tools.norm2(
            num.matrixmultiply(s, vector) 
            - value * num.matrixmultiply(m, vector))
