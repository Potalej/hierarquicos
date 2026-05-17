MODULE tree_mod
USE node_mod
IMPLICIT NONE
PRIVATE
PUBLIC TreeType

TYPE :: TreeType
    ! maximum depth (recursion)
    INTEGER :: max_depth = 32
    ! amplificator of the size of the quadtree-root
    REAL(8) :: side_amplificator = 1.2d0

    ! to save_txt
    INTEGER :: save_txt

    ! masses, positions and number of bodies (N)
    REAL(8), ALLOCATABLE :: m(:), x(:), y(:)
    INTEGER :: N

    ! root of the tree
    TYPE(NodeType), POINTER :: root => null()
CONTAINS
    PROCEDURE :: init
    PROCEDURE :: allocate_subnode, index_subnode, add_to_subnode, add
    PROCEDURE :: forces => evaluate_forces_over_p
    PROCEDURE :: evaluate_forces_of_node_over_p
END TYPE

CONTAINS

SUBROUTINE init (self, m, x, y, save_txt)
    CLASS(TreeType), INTENT(INOUT) :: self
    REAL(8), INTENT(IN) :: m(:), x(:), y(:)
    INTEGER, OPTIONAL :: save_txt
    REAL(8) :: infos_root(3)
    INTEGER :: p

    ! saving particles information
    self % N = SIZE(m)
    ALLOCATE(self % m(self % N))
    ALLOCATE(self % x(self % N))
    ALLOCATE(self % y(self % N))
    self % m = m
    self % x = x
    self % y = y

    ! init the root
    infos_root = node_size_center(x, y)
    ALLOCATE(self % root)
    CALL self % root % init(infos_root(1), infos_root(2), 1.2*infos_root(3), 0)

    ! if wants to save_txt
    self % save_txt = -1
    IF (PRESENT(save_txt)) THEN
        self % save_txt = save_txt
        WRITE (self % save_txt, *) self%root%cx, self%root%cy, self%root%halfside
    ENDIF

    ! add the children
    DO p = 1, self % N
        CALL self % add(self % root, p)
    END DO
END SUBROUTINE

SUBROUTINE allocate_subnode (self, node, index)
    CLASS(TreeType), INTENT(INOUT) :: self
    CLASS(NodeType),   INTENT(INOUT) :: node
    INTEGER,       INTENT(IN)    :: index
    REAL(8) :: h, h_half, cx, cy, cx_sub, cy_sub

    h = node % halfside
    h_half = h / 2.0d0
    cx = node % cx
    cy = node % cy

    IF (index == 0) THEN
        ALLOCATE(node % child_NO)
        cx_sub = cx - h_half
        cy_sub = cy + h_half
        CALL node % child_NO % init(cx_sub, cy_sub, h, node % depth + 1)
    ELSEIF (index == 1) THEN
        ALLOCATE(node % child_NE)
        cx_sub = cx + h_half
        cy_sub = cy + h_half
        CALL node % child_NE % init(cx_sub, cy_sub, h, node % depth + 1)
    ELSEIF (index == 2) THEN
        ALLOCATE(node % child_SO)
        cx_sub = cx - h_half
        cy_sub = cy - h_half
        CALL node % child_SO % init(cx_sub, cy_sub, h, node % depth + 1)
    ELSE
        ALLOCATE(node % child_SE)
        cx_sub = cx + h_half
        cy_sub = cy - h_half
        CALL node % child_SE % init(cx_sub, cy_sub, h, node % depth + 1)
    ENDIF

    IF (self % save_txt .NE. -1) WRITE (self % save_txt, *) node % depth + 1, cx_sub, cy_sub
END SUBROUTINE

SUBROUTINE index_subnode (self, node, p, index)
    CLASS(TreeType), INTENT(INOUT) :: self
    CLASS(NodeType),   INTENT(INOUT) :: node
    INTEGER,       INTENT(IN)    :: p
    INTEGER,       INTENT(INOUT) :: index
    LOGICAL :: top, right

    top   = (self % y(p) > node % cy)
    right = (self % x(p) > node % cx)

    IF (top) THEN
        IF (right) THEN
            index = 1
        ELSE
            index = 0
        ENDIF
    ELSE
        IF (right) THEN
            index = 3
        ELSE
            index = 2
        ENDIF
    ENDIF
    
END SUBROUTINE

SUBROUTINE add_to_subnode (self, node, p)
    CLASS(TreeType), INTENT(INOUT) :: self
    CLASS(NodeType), TARGET, INTENT(INOUT) :: node
    INTEGER, INTENT(IN) :: p
    CLASS(NodeType), POINTER :: subnode
    INTEGER :: subnode_index

    CALL self % index_subnode(node, p, subnode_index)
    CALL node % get_child_by_index(subnode_index, subnode)

    ! if the subnode isnt associated, so associate
    IF (.NOT. ASSOCIATED(subnode)) THEN
        CALL self % allocate_subnode(node, subnode_index)
    ENDIF

    CALL node % get_child_by_index(subnode_index, subnode)

    ! now add
    CALL self % add(subnode, p)
END SUBROUTINE

SUBROUTINE add (self, node, p)
    CLASS(TreeType), INTENT(INOUT) :: self
    CLASS(NodeType), INTENT(INOUT) :: node
    INTEGER,       INTENT(IN)    :: p ! particle index
    INTEGER :: old_p
    REAL(8) :: pm, px, py, old_mass

    ! get particle information
    pm = self % m(p)
    px = self % x(p)
    py = self % y(p)

    ! an empty node become a particle
    IF (node % particle == -1 .AND. node % is_leaf) THEN
        node % particle = p

        ! add directly the mass and center of mass
        node % mass = pm
        node % qcm_x = px
        node % qcm_y = py

        RETURN
    ENDIF

    ! update the mass and center of mass if the node isnt empty
    old_mass = node % mass
    node % mass = node % mass + pm
    node % qcm_x = (node % qcm_x * old_mass + px * pm) / node % mass
    node % qcm_y = (node % qcm_y * old_mass + py * pm) / node % mass

    ! if isnt empty, it become a twig
    IF (node % is_leaf) THEN
        ! we cannot go beyond the depth limit
        IF (node % depth >= self % max_depth) THEN
            PRINT *, "PROBLEMA SERIO !!!"
            STOP 0
        ENDIF

        ! in this case, its now a twig
        node % is_leaf = .FALSE.

        ! add the old particle as a particle per si
        old_p = node % particle
        node % particle = -1
        CALL self % add_to_subnode(node, old_p)
    ENDIF

    ! now add the new particle
    CALL self % add_to_subnode(node, p)
END SUBROUTINE

FUNCTION evaluate_forces_over_p (self, p, par_eps, par_theta, par_G) RESULT (forces)
    CLASS(TreeType), INTENT(INOUT) :: self
    INTEGER, INTENT(IN) :: p ! particle index
    REAL(8), INTENT(IN), OPTIONAL :: par_eps, par_theta, par_G ! parameters
    REAL(8) :: eps, theta, G, forces(2)

    ! default values
    eps = 0.0d0
    theta = 0.0d0
    G = 1.0d0

    ! replace if present
    IF (PRESENT(par_eps))   eps = par_eps
    IF (PRESENT(par_theta)) theta = par_theta
    IF (PRESENT(par_G))     G = par_G

    ! evaluate the forces over p
    forces = self % evaluate_forces_of_node_over_p(self % root, p, eps, theta, G)
END FUNCTION

RECURSIVE FUNCTION evaluate_forces_of_node_over_p (self, node, p, eps, theta, G) RESULT (forces)
    CLASS(TreeType), INTENT(INOUT) :: self
    CLASS(NodeType), INTENT(IN) :: node
    INTEGER, INTENT(IN) :: p ! particle index
    REAL(8), INTENT(IN) :: eps, theta, G
    REAL(8) :: forces(2) ! output
    REAL(8) :: dx, dy, dist2, dist, L, rab, f
    INTEGER :: i

    ! if no mass, no forces (it should not happen btw)
    IF (node % mass == 0) THEN
        forces = 0.0d0
        RETURN
    ENDIF

    ! if its the particle, no forces too
    IF (node % is_leaf .AND. node % particle == p) THEN
        forces = 0.0d0
        RETURN
    ENDIF

    ! in other cases, verify Barnes-Hut condition
    dx = node % qcm_x - self % x(p)
    dy = node % qcm_y - self % y(p)

    dist2 = dx*dx + dy*dy
    dist = SQRT(dist2)
    L = 2.0d0 * node % halfside

    ! if its a leaf or the BH condition is True
    IF (node % is_leaf .OR. L / dist < theta) THEN
        rab = SQRT(dist2 + eps*eps)
        f = G * self % m(p) * node % mass / (rab ** 3.0d0)
        forces(1) = f * dx
        forces(2) = f * dy
        RETURN
    ENDIF

    ! if its not a leaf nor the BH condition is valid, go to the leafs
    forces = 0.0d0
    IF (ASSOCIATED(node % child_NO)) THEN
        forces = forces + self % evaluate_forces_of_node_over_p(node % child_NO, p, eps, theta, G)
    ENDIF
    IF (ASSOCIATED(node % child_NE)) THEN
        forces = forces + self % evaluate_forces_of_node_over_p(node % child_NE, p, eps, theta, G)
    ENDIF
    IF (ASSOCIATED(node % child_SO)) THEN
        forces = forces + self % evaluate_forces_of_node_over_p(node % child_SO, p, eps, theta, G)
    ENDIF
    IF (ASSOCIATED(node % child_SE)) THEN
        forces = forces + self % evaluate_forces_of_node_over_p(node % child_SE, p, eps, theta, G)
    ENDIF
END FUNCTION

END MODULE