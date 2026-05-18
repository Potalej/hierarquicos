! api for compile with f2py
MODULE api_mod
    USE tree_mod
    IMPLICIT NONE

    PUBLIC
CONTAINS

SUBROUTINE generate_tree (m, x, y)
    REAL(8), INTENT(IN) :: m(:), x(:), y(:)
    TYPE(TreeType) :: tree

    CALL tree % init(m, x, y)
END SUBROUTINE

FUNCTION compute_forces (m, x, y, theta, G, eps) RESULT (forces)
    REAL(8), INTENT(IN) :: m(:), x(:), y(:)
    REAL(8), INTENT(IN) :: theta, G, eps
    
    TYPE(TreeType) :: tree
    INTEGER :: p
    REAL(8) :: forces(SIZE(m), 2)

    CALL tree % init(m, x, y)

    DO p = 1, SIZE(m)
        forces(p,:) = tree % forces(p, theta, G, eps)
    END DO
END FUNCTION

FUNCTION compute_forces_direct (m, x, y, G, eps) RESULT (forces)
    REAL(8), INTENT(IN) :: m(:), x(:), y(:)
    REAL(8), INTENT(IN) :: G, eps
    INTEGER :: a, b
    REAL(8) :: dx, dy, dist2, dist, f
    REAL(8) :: forces(SIZE(m), 2)

    forces = 0.0d0

    DO a = 2, SIZE(m)
        DO b = 1, a - 1
            dx = x(b) - x(a)
            dy = y(b) - y(a)
            dist2 = dx*dx + dy*dy + eps*eps
            dist = SQRT(dist2)
            f = G * m(a) * m(b) / (dist * dist2)

            ! forces over 'a'
            forces(a,1) = forces(a,1) + f * dx
            forces(a,2) = forces(a,2) + f * dy

            ! forces over 'b'
            forces(b,1) = forces(b,1) - f * dx
            forces(b,2) = forces(b,2) - f * dy
        END DO
    END DO
END FUNCTION

END MODULE