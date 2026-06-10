! api for compile with f2py (3d)
MODULE api_mod
    USE octree_mod
    IMPLICIT NONE

    PUBLIC
CONTAINS

SUBROUTINE generate_tree (m, x, y, z, quad)
    REAL(16), INTENT(IN) :: m(:), x(:), y(:), z(:)
    LOGICAL, OPTIONAL :: quad
    TYPE(OctreeType) :: tree

    CALL tree % init(m, x, y, z)
    IF (PRESENT(quad)) THEN
        IF (quad) CALL tree % evaluate_quad()
    ENDIF
END SUBROUTINE

FUNCTION compute_forces (m, x, y, z, theta, G, eps, quad) RESULT (forces)
    REAL(16), INTENT(IN) :: m(:), x(:), y(:), z(:)
    REAL(16), INTENT(IN) :: theta, G, eps
    LOGICAL, OPTIONAL :: quad
    
    TYPE(OctreeType) :: tree
    INTEGER :: p
    REAL(16) :: forces(SIZE(m), 3)
    REAL(16) :: theta2

    theta2 = theta * theta

    CALL tree % init(m, x, y, z)
    IF (PRESENT(quad)) THEN
        IF (quad) CALL tree % evaluate_quad()
    ENDIF

    !$OMP PARALLEL SHARED(forces) PRIVATE(p) 
    !$OMP DO
    DO p = 1, SIZE(m)
        forces(p,:) = tree % forces(p, theta2, G, eps)
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
END FUNCTION

FUNCTION compute_forces_direct (m, x, y, z, G, eps) RESULT (forces)
    REAL(16), INTENT(IN) :: m(:), x(:), y(:), z(:)
    REAL(16), INTENT(IN) :: G, eps
    INTEGER :: a, b
    REAL(16) :: dx, dy, dz, dist2, dist, f
    REAL(16) :: forces(SIZE(m), 3)

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

FUNCTION compute_forces_direct_par (m, x, y, z, G, eps) RESULT (forces)
    REAL(16), INTENT(IN) :: m(:), x(:), y(:), z(:)
    REAL(16), INTENT(IN) :: G, eps
    INTEGER :: a, b
    REAL(16) :: dx, dy, dz, dist2, dist, f
    REAL(16) :: forces(SIZE(m), 3)

    forces = 0.0d0

    !$OMP PARALLEL SHARED(forces) PRIVATE(a, b, dx, dy, dz, dist2, dist, f)
    !$OMP DO
    DO a = 1, SIZE(m)
        DO b = 1, SIZE(m)
            IF (a == b) CYCLE
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
        END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
END FUNCTION

FUNCTION total_energy (m, qs, ps, G, eps)
    REAL(16), INTENT(IN) :: m(:), qs(:,:), ps(:,:), G, eps
    REAL(16) :: total_energy, V, rab
    INTEGER :: a, b

    total_energy = 0.0_16
    DO a = 1, SIZE(m)
        total_energy = total_energy + 0.5_16*DOT_PRODUCT(ps(a,:),ps(a,:))/m(a)
        DO b = 1, a - 1
            rab = SQRT(NORM2(qs(b,:)-qs(a,:))**2 + eps*eps)
            V = - G * m(a) * m(b) / rab
            total_energy = total_energy + V
        END DO
    END DO 
END FUNCTION

FUNCTION total_linear_momentum (ps) RESULT(P)
    REAL(16), INTENT(IN) :: ps(:,:)
    REAL(16) :: P(3)
    INTEGER :: a

    P = 0.0_16
    DO a = 1, SIZE(ps,1)
        P = P + ps(a,:)
    END DO
END FUNCTION

FUNCTION cross_product (u, v)
  REAL(16), DIMENSION(:), INTENT(IN) :: u, v
  REAL(16), DIMENSION(3)             :: cross_product

  cross_product(1) =  u(2)*v(3)-v(2)*u(3)
  cross_product(2) = -u(1)*v(3)+v(1)*u(3)
  cross_product(3) =  u(1)*v(2)-v(1)*u(2)
END FUNCTION

FUNCTION total_angular_momentum (qs, ps) RESULT(J)
    REAL(16), INTENT(IN) :: qs(:,:), ps(:,:)
    REAL(16) :: J(3)
    INTEGER :: a

    J = 0.0_16
    DO a = 1, SIZE(ps,1)
        J = J + cross_product(qs(a,:), ps(a,:))
    END DO
END FUNCTION

END MODULE