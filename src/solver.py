import math

import pylinear.array as num
import pylinear.operation as op

import stopwatch
import element
import mesh_function





def solveSPDSystem(matrix_op, rhs, start_vector = None):
    h,w = matrix_op.shape
    if start_vector is None:
        x = num.zeros((w,), matrix_op.typecode())
    else:
        x = start_vector.copy()

    job = stopwatch.Job("solve")
    matrix_inv_op = op.CGOperator.make(matrix_op)
    matrix_inv_op.apply(rhs, x)
    job.done()
    print "  iterations:", matrix_inv_op.last_iteration_count

    residual = x.copy()
    matrix_op.apply(x, residual)
    residual -= rhs

    print "  absolute residual: ", op.norm_2(residual)
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

    job = stopwatch.Job("stiffness matrix")
    for el in mesh.elements():
        nonus = [number_assignment[node] for node in el.nodes()]
        s.add_scattered(
            nonus, nonus,
            num.asarray(el.getVolumeIntegralsOverDifferentiatedFormFunctions(),
                        typecode))

    job.done()

    return s




def buildMassMatrix(mesh, number_assignment, weight_function, typecode = num.Float):
    dof_count = len(number_assignment)
    m = num.zeros((dof_count, dof_count), typecode, num.SparseBuildMatrix)

    job = stopwatch.Job("mass matrix")
    for el in mesh.elements():
        nonus = [number_assignment[node] for node in el.nodes()]
        m.add_scattered(
            nonus, nonus,
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

    job = stopwatch.Job("rhs")
    for el in mesh.elements():
        this_el_b = el.getVolumeIntegralsOverFormFunction(
            lambda x,formfunc_value: f(x) * formfunc_value)
        for node, v in zip(el.nodes(), this_el_b):
            try:
                b[number_assignment[node]] += v
            except KeyError:
                pass
    job.done()

    job = stopwatch.Job("bcs: dirichlet")
    for node, (boundary_value, other_nodes) in node_constraints.iteritems():
        assert len(other_nodes) == 0
        if boundary_value != 0:
            nonu = complete_number_assignment[node]
            b += boundary_value * full_s[:dof_count, nonu]
    job.done()

    #visualization.writeGnuplotSparsityPattern(",,s.data", s)

    compiled_s = num.asarray(s, s.typecode(), num.SparseExecuteMatrix)
    s_op = op.MatrixOperator.make(compiled_s)
    
    complete_vec = mesh_function.makeCompleteVector(
        complete_number_assignment,
        solveSPDSystem(s_op, -b),
        node_constraints)
    return mesh_function.tMeshFunction(
        mesh, complete_number_assignment, complete_vec)




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
    return (vec*s*vec.H)/(vec*m*vec.H)




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
                                
        self.FullS = buildStiffnessMatrix(mesh, complete_number_assignment, typecode)
        self.FullM = buildMassMatrix(mesh, complete_number_assignment, g, typecode)

        if f is not None:
            job = stopwatch.Job("f")
    
            for el in mesh.elements():
                nonus = [number_assignment[node] for node in el.nodes()]
                self.FullS.add_scattered(
                    nonus, nonus,
                    num.asarray(el.getVolumeIntegralsOverFormFunctions(f), typecode))
            job.done()

        self.FullSEx = num.asarray(self.FullS, flavor=num.SparseExecuteMatrix)
        self.FullSOp = op.MatrixOperator.make(self.FullSEx)

        self.FullMEx = num.asarray(self.FullM, flavor=num.SparseExecuteMatrix)
        self.FullMOp = op.MatrixOperator.make(self.FullMEx)

        self.Constraints = None
        self.ConstrainedSEx = self.ConstrainedMEx = None
        self.ConstraintOp = self.ConstraintHermOp = None

    def stiffnessMatrix(self):
        return self.FullSEx

    def massMatrix(self):
        return self.FullMEx

    def nodeNumberAssignment(self):
        return self.CompleteNumberAssignment

    def currentConstraints(self):
        return self.Constraints

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
        
        a_ex = num.asarray(a, flavor=num.SparseExecuteMatrix)

        job = stopwatch.Job("constrained matrices")
        self.ConstrainedSEx = a_ex * self.FullSEx * a_ex.H
        self.ConstrainedMEx = a_ex * self.FullMEx * a_ex.H
        job.done()
        self.Constraints = constraints

        return a

    def solve(self, sigma, number_of_eigenvalues = 20, tolerance = 1e-10,
              warning_threshold = 10):
        if self.Constraints is None:
            raise RuntimeError, "need to set up constraints before solving"

        m_op = op.MatrixOperator.make(self.ConstrainedMEx)

        if sigma == 0:
            m_inv_op = op.CGOperator.make(
                m_op, 
                tolerance=tolerance,
                precon_op=op.SSORPreconditioner.make(self.ConstrainedMEx))

            job = stopwatch.Job("arpack rci")
            eigen_result = op.operator_eigenvectors(
                m_inv_op * op.MatrixOperator.make(self.ConstrainedSEx), 
                number_of_eigenvalues,
                m_op, 
                tolerance=tolerance,
                which_eigenvalues=op.SMALLEST_MAGNITUDE)
            job.done()
        else:
            shifted_mat_ex = self.ConstrainedSEx - sigma*self.ConstrainedMEx
            job = stopwatch.Job("arpack rci")
            eigen_result = op.operator_eigenvectors(
                op.UMFPACKOperator.make(shifted_mat_ex) * m_op, 
                number_of_eigenvalues,
                m_op, 
                sigma, 
                tolerance=tolerance)
            job.done()
  
        result = [(eigenvalue, mesh_function.tMeshFunction(
            self.Mesh,
            self.CompleteNumberAssignment,
            mesh_function.makeCompleteVector(self.CompleteNumberAssignment,
                                             num.conjugate(eigenvector),
                                             self.Constraints)))
                for eigenvalue, eigenvector in eigen_result]

        # perform some invariant checking
        if warning_threshold is not None:
            for evalue, mf in result:
                if self.computeEigenpairResidual(evalue, mf) > warning_threshold * tolerance:
                    print "PRECSION WARNING for ARPACK output."

        return result

    def computeEigenpairResidual(self, eigenvalue, mesh_func):
        if self.Constraints is None:
            raise RuntimeError, "need to set up constraints before verifying"

        dof_count = len(self.Mesh.dofManager())- len(self.Constraints)
        vector = num.conjugate(mesh_func.vector()[:dof_count])

        def m_sp(x,y): return x * (self.ConstrainedMEx*y.H)
        def m_norm(x): return math.sqrt(abs(m_sp(x,x)))

        return m_norm(
            self.ConstrainedSEx*vector - eigenvalue * self.ConstrainedMEx*vector)
