MODULE node_mod
IMPLICIT NONE
PRIVATE
PUBLIC NodeType, node_size_center

! Node
TYPE :: NodeType
    ! geometric center
    REAL(8) :: cx, cy
    ! square side halfsize
    REAL(8) :: halfside

    ! mass and center of mass
    REAL(8) :: mass, qcm_x, qcm_y

    ! particle index
    INTEGER :: particle

    ! depth in the tree
    INTEGER :: depth

    ! children
    TYPE(NodeType), POINTER :: child_NO => null()
    TYPE(NodeType), POINTER :: child_NE => null()
    TYPE(NodeType), POINTER :: child_SO => null()
    TYPE(NodeType), POINTER :: child_SE => null()

    ! if its a leaf
    LOGICAL :: is_leaf
CONTAINS
    PROCEDURE :: init, get_child_by_index
END TYPE

CONTAINS

SUBROUTINE init (self, cx, cy, side, depth)
    CLASS(NodeType), INTENT(INOUT) :: self
    REAL(8), INTENT(IN) :: cx, cy, side
    INTEGER, INTENT(IN), OPTIONAL :: depth

    self % cx = cx
    self % cy = cy

    self % halfside = side / 2.0d0

    self % mass  = 0.0d0
    self % qcm_x = 0.0d0
    self % qcm_y = 0.0d0

    self % particle = -1
    IF (PRESENT(depth)) THEN
        self % depth = depth
    ELSE
        self % depth = 0
    ENDIF

    self % is_leaf = .TRUE.
END SUBROUTINE

SUBROUTINE get_child_by_index (self, index, subnode)
    CLASS(NodeType), TARGET, INTENT(INOUT) :: self
    CLASS(NodeType), POINTER, INTENT(OUT) :: subnode
    INTEGER,   INTENT(IN)    :: index

    IF (index == 0) THEN
        subnode => self % child_NO
    ELSEIF (index == 1) THEN
        subnode => self % child_NE
    ELSEIF (index == 2) THEN
        subnode => self % child_SO
    ELSE
        subnode => self % child_SE
    ENDIF
END SUBROUTINE

FUNCTION node_size_center (x, y) RESULT (infos)
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

END MODULE