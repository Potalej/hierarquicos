MODULE forces_mod
    IMPLICIT NONE
    PRIVATE
    PUBLIC compute_forces_direct
    INTEGER, PARAMETER  :: pf  = SELECTED_REAL_KIND(15, 307)
CONTAINS

FUNCTION compute_forces_direct (m, x, y, z, G, eps2, nt) RESULT (forces)
    REAL(pf), INTENT(IN) :: m(:), x(:), y(:), z(:)
    REAL(pf), INTENT(IN) :: G, eps2
    INTEGER, INTENT(IN) :: nt
    INTEGER :: a, b
    REAL(pf) :: dx, dy, dz, dist2, dist, f
    REAL(pf) :: forces(3, SIZE(m))

    forces = 0.0_pf

    IF (nt == 1) THEN
        DO a = 2, SIZE(m)
            DO b = 1, a - 1
                dx = x(b) - x(a)
                dy = y(b) - y(a)
                dz = z(b) - z(a)
                dist2 = dx*dx + dy*dy + dz*dz + eps2
                dist = SQRT(dist2)
                f = G * m(a) * m(b) / (dist * dist2)

                ! forces over 'a'
                forces(1,a) = forces(1,a) + f * dx
                forces(2,a) = forces(2,a) + f * dy
                forces(3,a) = forces(3,a) + f * dz

                ! forces over 'b'
                forces(1,b) = forces(1,b) - f * dx
                forces(2,b) = forces(2,b) - f * dy
                forces(3,b) = forces(3,b) - f * dz
            END DO
        END DO
    ELSE
        !$OMP PARALLEL DO PRIVATE(a,b,dx,dy,dz,dist2,dist,f) NUM_THREADS(nt) &
        !$OMP SCHEDULE(DYNAMIC) REDUCTION(+:forces)
        DO a = 2, SIZE(m)
            DO b = 1, a - 1
                dx = x(b) - x(a)
                dy = y(b) - y(a)
                dz = z(b) - z(a)
                dist2 = dx*dx + dy*dy + dz*dz + eps2
                dist = SQRT(dist2)
                f = G * m(a) * m(b) / (dist * dist2)

                ! forces over 'a'
                forces(1,a) = forces(1,a) + f * dx
                forces(2,a) = forces(2,a) + f * dy
                forces(3,a) = forces(3,a) + f * dz

                ! forces over 'b'
                forces(1,b) = forces(1,b) - f * dx
                forces(2,b) = forces(2,b) - f * dy
                forces(3,b) = forces(3,b) - f * dz
            END DO
        END DO
        !$OMP END PARALLEL DO
    ENDIF
END FUNCTION

END MODULE