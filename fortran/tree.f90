MODULE tree_mod
IMPLICIT NONE
PRIVATE
PUBLIC TreeType

TYPE :: TreeType
    ! max depth of the tree
    INTEGER :: max_depth = 32
    ! amplificator of the size of the quadtree-root
    REAL(8) :: side_amplificator = 1.2d0

    ! to save_txt
    INTEGER :: save_txt

    ! masses, positions and number of bodies (N)
    REAL(8), ALLOCATABLE :: m(:), x(:), y(:)
    INTEGER :: N

    ! nodes
    INTEGER :: number_of_nodes = 0
    INTEGER :: max_number_of_nodes
    REAL(8), ALLOCATABLE :: ns_cx(:), ns_cy(:) ! centers
    REAL(8), ALLOCATABLE :: ns_halfside(:), ns_L2(:)
    REAL(8), ALLOCATABLE :: ns_mass(:), ns_qcm_x(:), ns_qcm_y(:)
    INTEGER, ALLOCATABLE :: ns_particle(:)
    LOGICAL, ALLOCATABLE :: ns_is_leaf(:)
    INTEGER, ALLOCATABLE :: ns_child(:,:)
    INTEGER, ALLOCATABLE :: ns_depth(:)
CONTAINS
    PROCEDURE :: init
    PROCEDURE :: allocate_nodes, add_node, allocate_subnode, index_subnode, add_to_subnode, add
    PROCEDURE :: forces => evaluate_forces_over_p
END TYPE

CONTAINS

SUBROUTINE init (self, m, x, y, save_txt)
! this subroutine inits the tree by allocating the global vectors and adding each particle
! in a node. if its the case it saves the root information too.
    CLASS(TreeType), INTENT(INOUT) :: self
    REAL(8), INTENT(IN) :: m(:), x(:), y(:)
    INTEGER, OPTIONAL :: save_txt
    REAL(8) :: infos_root(3)
    INTEGER :: p, idx_root

    ! saving particles information
    self % N = SIZE(m)
    ALLOCATE(self % m(self % N))
    ALLOCATE(self % x(self % N))
    ALLOCATE(self % y(self % N))
    self % m = m
    self % x = x
    self % y = y

    ! init the root
    self % max_number_of_nodes = 8 * self % N
    self % number_of_nodes = 0
    CALL self % allocate_nodes()
    
    infos_root = node_size_center(x, y)
    CALL self % add_node(infos_root(1), infos_root(2), 1.2*infos_root(3), 0, idx_root)

    ! if wants to save_txt
    self % save_txt = -1
    IF (PRESENT(save_txt)) THEN
        self % save_txt = save_txt
        WRITE (self % save_txt, *) self%ns_cx(1), self%ns_cy(1), self%ns_halfside(1)
    ENDIF

    ! add the children to root
    DO p = 1, self % N
        CALL self % add(idx_root, p)
    END DO
END SUBROUTINE

SUBROUTINE allocate_nodes (self)
! this subroutine allocates the global vectors that stores informations about the nodes.
! if the tree already have nodes, so its a reallocation. in that case, the old information
! are stored in temp vectors and the vectors are reallocated with twice times the old size
! without lose the old information.
! @calledby init, add_node
    CLASS(TreeType), INTENT(INOUT) :: self
    REAL(8), ALLOCATABLE :: temp_real(:)
    INTEGER, ALLOCATABLE :: temp_int(:), temp_int_2(:,:)
    LOGICAL, ALLOCATABLE :: temp_logical(:)
    INTEGER :: old_size, new_size

    ! if the tree already exists, so is the case of reallocation
    IF (self % number_of_nodes > 0) THEN
        old_size = self % max_number_of_nodes
        self % max_number_of_nodes = 2 * old_size
        new_size = self % max_number_of_nodes
        
        ! allocate the temp vectors
        ALLOCATE(temp_real(old_size))
        ALLOCATE(temp_int(old_size))
        ALLOCATE(temp_int_2(4,old_size))
        ALLOCATE(temp_logical(old_size))

        ! now deallocate and reallocate
        temp_real = self % ns_cx
        DEALLOCATE(self % ns_cx)
        ALLOCATE(self % ns_cx(new_size))
        self % ns_cx(1:new_size) = temp_real

        temp_real = self % ns_cy
        DEALLOCATE(self % ns_cy)
        ALLOCATE(self % ns_cy(new_size))
        self % ns_cy(1:new_size) = temp_real

        temp_real = self % ns_halfside
        DEALLOCATE(self % ns_halfside)
        ALLOCATE(self % ns_halfside(new_size))
        self % ns_halfside(1:new_size) = temp_real

        temp_real = self % ns_L2
        DEALLOCATE(self % ns_L2)
        ALLOCATE(self % ns_L2(new_size))
        self % ns_L2(1:new_size) = temp_real

        temp_real = self % ns_mass
        DEALLOCATE(self % ns_mass)
        ALLOCATE(self % ns_mass(new_size))
        self % ns_mass(1:new_size) = temp_real

        temp_real = self % ns_qcm_x
        DEALLOCATE(self % ns_qcm_x)
        ALLOCATE(self % ns_qcm_x(new_size))
        self % ns_qcm_x(1:new_size) = temp_real

        temp_real = self % ns_qcm_y
        DEALLOCATE(self % ns_qcm_y)
        ALLOCATE(self % ns_qcm_y(new_size))
        self % ns_qcm_y(1:new_size) = temp_real

        temp_int = self % ns_particle
        DEALLOCATE(self % ns_particle)
        ALLOCATE(self % ns_particle(new_size))
        self % ns_particle(1:new_size) = temp_int

        temp_logical = self % ns_is_leaf
        DEALLOCATE(self % ns_is_leaf)
        ALLOCATE(self % ns_is_leaf(new_size))
        self % ns_is_leaf(1:new_size) = temp_logical
        
        temp_int_2 = self % ns_child
        DEALLOCATE(self % ns_child)
        ALLOCATE(self % ns_child(4,new_size))
        self % ns_child(:,1:new_size) = temp_int_2

        DEALLOCATE(temp_real, temp_int, temp_int_2, temp_logical)
    ELSE
        ALLOCATE(self % ns_cx(self % max_number_of_nodes))
        ALLOCATE(self % ns_cy(self % max_number_of_nodes))
        ALLOCATE(self % ns_halfside(self % max_number_of_nodes))
        ALLOCATE(self % ns_L2(self % max_number_of_nodes))
        ALLOCATE(self % ns_mass(self % max_number_of_nodes))
        ALLOCATE(self % ns_qcm_x(self % max_number_of_nodes))
        ALLOCATE(self % ns_qcm_y(self % max_number_of_nodes))
        ALLOCATE(self % ns_particle(self % max_number_of_nodes))
        ALLOCATE(self % ns_is_leaf(self % max_number_of_nodes))
        ALLOCATE(self % ns_depth(self % max_number_of_nodes))
        ALLOCATE(self % ns_child(4, self % max_number_of_nodes))
    ENDIF
END SUBROUTINE

FUNCTION node_size_center (x, y) RESULT (infos)
! this subroutine gets the node size and center (in plane xy) based in the position
! vectors, giving the minimal square that contains every body.
! @calledby init
    REAL(8), INTENT(IN) :: x(:), y(:)
    REAL(8) :: xmin, xmax, ymin, ymax
    REAL(8) :: infos(3)

    xmin = MINVAL(x)
    xmax = MAXVAL(x)
    ymin = MINVAL(y)
    ymax = MAXVAL(y)

    infos(1) = 0.5d0 * (xmin + xmax)
    infos(2) = 0.5d0 * (ymin + ymax)
    infos(3) = MAXVAL((/ xmax - xmin, ymax - ymin /))
END FUNCTION

SUBROUTINE add_node (self, cx, cy, side, depth, idx)
! this subroutine adds a node to the tree as a leaf. if its necessary, it reallocate the
! global node state vectors with twice the original size.
! @calledby init, allocate_subnode
    CLASS(TreeType), INTENT(INOUT) :: self
    REAL(8), INTENT(IN) :: cx, cy, side
    INTEGER, INTENT(IN) :: depth
    INTEGER, INTENT(INOUT) :: idx

    self % number_of_nodes = self % number_of_nodes + 1
    idx = self % number_of_nodes
    
    IF (idx > self % max_number_of_nodes) THEN
        CALL self % allocate_nodes()
        PRINT *, '[DEBUG] REALLOCATING!'
    END IF

    ! starts without children and being a leaf
    self % ns_child(:, idx) = -1
    self % ns_is_leaf(idx) = .TRUE.
    self % ns_depth(idx) = depth

    self % ns_cx(idx) = cx
    self % ns_cy(idx) = cy
    self % ns_halfside(idx) = side / 2.0d0
    self % ns_L2(idx) = side * side
    self % ns_mass(idx) = 0.0d0
    self % ns_qcm_x(idx) = 0.0d0
    self % ns_qcm_y(idx) = 0.0d0
    self % ns_particle(idx) = -1
END SUBROUTINE

SUBROUTINE allocate_subnode (self, node_idx, index)
! this subroutine allocates a subnode based in the index of the subnode wrt the node.
! the index is the default: 1-NO, 2-NE, 3-SO, 4-SE.
! @calledby add_to_subnode
    CLASS(TreeType), INTENT(INOUT) :: self
    INTEGER, INTENT(IN) :: node_idx
    INTEGER, INTENT(IN) :: index
    INTEGER :: subnode_idx, d
    REAL(8) :: h, h_half, cx, cy, cx_sub, cy_sub
    
    h = self % ns_halfside(node_idx)
    h_half = h / 2.0d0
    cx = self % ns_cx(node_idx)
    cy = self % ns_cy(node_idx)
    d = self % ns_depth(node_idx) + 1

    ! NO
    IF (index == 1) THEN
        cx_sub = cx - h_half
        cy_sub = cy + h_half
        ! create subnode
        CALL self % add_node(cx_sub, cy_sub, h, d, subnode_idx)
        ! allocate it as child of node
        self % ns_child(1,node_idx) = subnode_idx
    ! NE
    ELSEIF (index == 2) THEN
        cx_sub = cx + h_half
        cy_sub = cy + h_half
        ! create subnode
        CALL self % add_node(cx_sub, cy_sub, h, d, subnode_idx)
        ! allocate it as child of node
        self % ns_child(2,node_idx) = subnode_idx
    ! SO
    ELSEIF (index == 3) THEN
        cx_sub = cx - h_half
        cy_sub = cy - h_half
        ! create subnode
        CALL self % add_node(cx_sub, cy_sub, h, d, subnode_idx)
        ! allocate it as child of node
        self % ns_child(3,node_idx) = subnode_idx
    ! SE
    ELSE
        cx_sub = cx + h_half
        cy_sub = cy - h_half
        ! create subnode
        CALL self % add_node(cx_sub, cy_sub, h, d, subnode_idx)
        ! allocate it as child of node
        self % ns_child(4,node_idx) = subnode_idx
    ENDIF

    IF (self % save_txt .NE. -1) WRITE (self % save_txt, *) d, cx_sub, cy_sub
END SUBROUTINE

SUBROUTINE index_subnode (self, node_idx, p, index)
! based in the position of the particle it returns a quadrant wrt to the node: 
! 1-NO, 2-NE, 3-SO, 4-SE.
! @calledby add_to_subnode
    CLASS(TreeType), INTENT(INOUT) :: self
    INTEGER, INTENT(IN) :: node_idx
    INTEGER, INTENT(IN)    :: p
    INTEGER, INTENT(INOUT) :: index
    LOGICAL :: top, right

    top   = (self % y(p) > self % ns_cy(node_idx))
    right = (self % x(p) > self % ns_cx(node_idx))

    IF (top) THEN
        IF (right) THEN
            index = 2
        ELSE
            index = 1
        ENDIF
    ELSE
        IF (right) THEN
            index = 4
        ELSE
            index = 3
        ENDIF
    ENDIF
    
END SUBROUTINE

SUBROUTINE add_to_subnode (self, node_idx, p)
! given a node/twig and a particle, this adds the particle to the correct leaf based
! in its position.
! @calledby add
    CLASS(TreeType), INTENT(INOUT) :: self
    INTEGER, INTENT(IN) :: node_idx
    INTEGER, INTENT(IN) :: p
    INTEGER :: subnode
    INTEGER :: subnode_index

    CALL self % index_subnode(node_idx, p, subnode_index)
    subnode = self % ns_child(subnode_index, node_idx)

    ! if the subnode isnt associated
    IF (subnode == -1) THEN
        CALL self % allocate_subnode(node_idx, subnode_index)
    ENDIF

    subnode_index = self % ns_child(subnode_index, node_idx)

    ! now add
    CALL self % add(subnode_index, p)
END SUBROUTINE

SUBROUTINE add (self, node_idx, p)
! given a node and a particle, it adds the particle to the node as a particle if the
! node is empty and as a leaf if the node is already a leaf/particle.
! @calledby init
    CLASS(TreeType), INTENT(INOUT) :: self
    INTEGER, INTENT(IN) :: node_idx
    INTEGER, INTENT(IN) :: p ! particle index
    INTEGER :: old_p
    REAL(8) :: pm, px, py, old_mass

    ! get particle information
    pm = self % m(p)
    px = self % x(p)
    py = self % y(p)

    ! an empty node become a particle
    IF (self % ns_particle(node_idx) == -1 .AND. self % ns_is_leaf(node_idx)) THEN
        self % ns_particle(node_idx) = p

        ! add directly the mass and center of mass
        self % ns_mass(node_idx) = pm
        self % ns_qcm_x(node_idx) = px
        self % ns_qcm_y(node_idx) = py

        RETURN
    ENDIF

    ! update the mass and center of mass if the node isnt empty
    old_mass = self % ns_mass(node_idx)
    self % ns_mass(node_idx) = self % ns_mass(node_idx) + pm
    self % ns_qcm_x(node_idx) = (self % ns_qcm_x(node_idx) * old_mass + px * pm) / self % ns_mass(node_idx)
    self % ns_qcm_y(node_idx) = (self % ns_qcm_y(node_idx) * old_mass + py * pm) / self % ns_mass(node_idx)

    ! if isnt empty, it become a twig
    IF (self % ns_is_leaf(node_idx)) THEN
        ! we cannot go beyond the depth limit
        IF (self % ns_depth(node_idx) >= self % max_depth) THEN
            PRINT *, "BIG PROBLEM !!! MAX DEPTH !!!"
            STOP 0
        ENDIF

        ! in this case, its now a twig
        self % ns_is_leaf(node_idx) = .FALSE.

        ! add the old particle as a particle per si
        old_p = self % ns_particle(node_idx)
        self % ns_particle(node_idx) = -1
        CALL self % add_to_subnode(node_idx, old_p)
    ENDIF

    ! now add the new particle
    CALL self % add_to_subnode(node_idx, p)
END SUBROUTINE

FUNCTION evaluate_forces_over_p (self, p, par_theta2, par_G, par_eps) RESULT (forces)
! given a particle and the parameters, this evaluates the forces over the particle
! using the Barnes-Hut criterion. it uses the depth first traversal (DFS) algorithm 
! to evaluate the forces in the tree.
    CLASS(TreeType), INTENT(INOUT) :: self
    INTEGER, INTENT(IN) :: p ! particle index
    REAL(8), INTENT(IN), OPTIONAL :: par_theta2, par_G, par_eps ! parameters
    REAL(8) :: eps, theta2, G

    REAL(8) :: pm, px, py
    REAL(8) :: forces(2)
    REAL(8) :: dx, dy, dist2, L2, f, rab

    INTEGER :: stack(self % number_of_nodes)
    INTEGER :: top, node_idx, i, child_idx

    ! default values
    eps = 0.0d0
    theta2 = 0.0d0
    G = 1.0d0

    ! replace if present
    IF (PRESENT(par_eps))   eps = par_eps
    IF (PRESENT(par_theta2)) theta2 = par_theta2
    IF (PRESENT(par_G))     G = par_G

    ! particle cache info
    pm = self % m(p)
    px = self % x(p)
    py = self % y(p)

    ! initialize
    forces = 0.0d0
    top = 1
    stack(top) = 1

    ! dfs iterative
    DO WHILE (top > 0)

        node_idx = stack(top)
        top = top - 1

        ! empty node
        IF (self % ns_mass(node_idx) == 0.0d0) CYCLE

        ! self interaction
        IF (self % ns_is_leaf(node_idx)) THEN
            IF (self % ns_particle(node_idx) == p) CYCLE
        ENDIF

        ! geometry
        dx = self % ns_qcm_x(node_idx) - px
        dy = self % ns_qcm_y(node_idx) - py
        dist2 = dx*dx + dy*dy
        L2 = self % ns_L2(node_idx)

        ! leaf or bh criterion
        IF (self % ns_is_leaf(node_idx) .OR. L2 < theta2 * dist2) THEN
            rab = SQRT(dist2 + eps*eps)
            f = G * pm * self % ns_mass(node_idx) / (rab * rab * rab)
            forces(1) = forces(1) + f * dx
            forces(2) = forces(2) + f * dy
        ! if isnt leaf nor bh is valid, so go to children
        ELSE
            ! push children in the vec
            DO i = 1, 4
                child_idx = self % ns_child(i, node_idx)
                IF (child_idx .NE. -1) THEN
                    top = top + 1
                    stack(top) = child_idx
                ENDIF
            END DO
        ENDIF
    END DO
END FUNCTION

END MODULE