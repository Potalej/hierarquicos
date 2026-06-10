! api for compile with f2py (3d)
MODULE api_mod
    USE octree_mod
    IMPLICIT NONE

    PUBLIC
CONTAINS

SUBROUTINE generate_tree (m, x, y, z, quad)
    REAL(8), INTENT(IN) :: m(:), x(:), y(:), z(:)
    LOGICAL, OPTIONAL :: quad
    TYPE(OctreeType) :: tree

    CALL tree % init(m, x, y, z)
    IF (PRESENT(quad)) THEN
        IF (quad) CALL tree % evaluate_quad()
    ENDIF
END SUBROUTINE

FUNCTION compute_forces (m, x, y, z, theta, G, eps, quad) RESULT (forces)
    REAL(8), INTENT(IN) :: m(:), x(:), y(:), z(:)
    REAL(8), INTENT(IN) :: theta, G, eps
    LOGICAL, OPTIONAL :: quad
    
    TYPE(OctreeType) :: tree
    INTEGER :: p
    REAL(8) :: forces(SIZE(m), 3)

    CALL tree % init(m, x, y, z)
    IF (PRESENT(quad)) THEN
        IF (quad) CALL tree % evaluate_quad()
    ENDIF

    !$OMP PARALLEL SHARED(forces) PRIVATE(p) 
    !$OMP DO
    DO p = 1, SIZE(m)
        forces(p,:) = tree % forces(p, theta, G, eps)
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
END FUNCTION

FUNCTION compute_forces_direct (m, x, y, z, G, eps) RESULT (forces)
    REAL(8), INTENT(IN) :: m(:), x(:), y(:), z(:)
    REAL(8), INTENT(IN) :: G, eps
    INTEGER :: a, b
    REAL(8) :: dx, dy, dz, dist2, dist, f
    REAL(8) :: forces(SIZE(m), 3)

    forces = 0.0d0

    DO a = 2, SIZE(m)
        DO b = 1, a - 1
            dx = x(b) - x(a)
            dy = y(b) - y(a)
            dz = z(b) - z(a)
            dist2 = dx*dx + dy*dy + dz*dz + eps*eps
            dist = SQRT(dist2)
            f = G * m(a) * m(b) / (dist * dist2)

            ! forces over 'a'
            forces(a,1) = forces(a,1) + f * dx
            forces(a,2) = forces(a,2) + f * dy
            forces(a,3) = forces(a,3) + f * dz

            ! forces over 'b'
            forces(b,1) = forces(b,1) - f * dx
            forces(b,2) = forces(b,2) - f * dy
            forces(b,3) = forces(b,3) - f * dz
        END DO
    END DO
END FUNCTION

END MODULE