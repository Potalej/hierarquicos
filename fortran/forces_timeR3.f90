! # TREE TIME
! 
! Program to test the time to construct the tree
! with N in [10, 1e4] (10 by 10) with 50 tests
! by value of N.
! It achieves the O(N log N), as expected :)
PROGRAM forces_time
    USE octree_mod
    USE omp_lib
    IMPLICIT NONE

    INTEGER, PARAMETER :: pf  = SELECTED_REAL_KIND(15, 307)

    INTEGER :: file, nt
    REAL(pf) :: thetas(5), eps2

    INTEGER :: Nmin, Nmax, Nstep
    LOGICAL :: quad

    file = 45
    nt = 10 ! number of threads
    ! eps2 = (0.1_pf)**2
    eps2 = 0.0_pf

    thetas(1) = 0.0_pf
    thetas(2) = 0.25_pf
    thetas(3) = 0.5_pf
    thetas(4) = 0.75_pf
    thetas(5) = 1.0_pf

    ! thetas(1) = 0.2_pf
    ! thetas(2) = 0.25_pf
    ! thetas(3) = 0.5_pf
    ! thetas(4) = 0.75_pf
    ! thetas(5) = 1.5_pf

    ! thetas(1) = 0.2_pf
    ! thetas(2) = 0.3_pf
    ! thetas(3) = 0.4_pf
    ! thetas(4) = 0.5_pf
    ! thetas(5) = 0.6_pf
    ! thetas(6) = 0.7_pf
    ! thetas(7) = 0.8_pf
    ! thetas(8) = 0.9_pf
    ! thetas(9) = 1.0_pf

    Nmin = 500
    Nmax = 10000
    Nstep = 500

    !!OPEN(file, file = "out/forces_time_morton_parallel_20C_N10K.txt", status="replace")
    ! OPEN(file, file = "out/test_morton_forces_time.txt", status="replace")
    ! CALL test_forces_time_morton(Nmin, Nmax, Nstep, thetas, eps2, 20, file, nt, quad)
    ! CLOSE(file)

    !!OPEN(file, file = "forces_time_without_morton_sequential.txt", status="replace")
    print *, '# monopole'
    OPEN(file, file = "out/test_multipoles_thetas_monopole.txt", status="replace")
    CALL test_forces_time(Nmin, Nmax, Nstep, thetas, eps2, 20, file, nt, 1)

    print *, '# quadrupole'
    OPEN(file, file = "out/test_multipoles_thetas_quadrupole.txt", status="replace")
    CALL test_forces_time(Nmin, Nmax, Nstep, thetas, eps2, 20, file, nt, 4)

    print *, '# octupole'
    OPEN(file, file = "out/test_multipoles_thetas_octupole.txt", status="replace")
    CALL test_forces_time(Nmin, Nmax, Nstep, thetas, eps2, 20, file, nt, 8)
    CLOSE(file)
CONTAINS

SUBROUTINE generate_initial_values (N, m, x, y, z)
    INTEGER, INTENT(IN) :: N
    REAL(pf), INTENT(INOUT) :: m(N), x(N), y(N), z(N)
    REAL(pf) :: qs(3,N)

    CALL RANDOM_SEED()
    m = 1.0_pf

    CALL RANDOM_NUMBER(qs)
    qs = 10.0_pf * (2.0_pf * qs - 1.0_pf)

    x = qs(1,:)
    y = qs(2,:)
    z = qs(3,:)
END SUBROUTINE

SUBROUTINE test_forces_time (Nmin, Nmax, Nstep, thetas, eps2, tests, file, nt, multipole)
    INTEGER,  INTENT(IN) :: Nmin, Nmax, Nstep, tests, file, nt
    REAL(pf), INTENT(IN) :: thetas(:), eps2
    INTEGER,  INTENT(IN) :: multipole
    INTEGER :: N, i_test, i_theta, p, timer
    REAL(pf), ALLOCATABLE :: m(:), x(:), y(:), z(:), forces(:,:), forces0(:,:)
    REAL(pf) :: time_start, time_finish, total, time_generate_tree
    REAL(pf) :: theta, erro, time_quad
    CLASS(OctreeType), ALLOCATABLE :: tree
    REAL(pf) :: a_scale

    DO N = Nmin, Nmax, Nstep
        PRINT *, 'N=', N
        ALLOCATE(m(N))
        ALLOCATE(x(N))
        ALLOCATE(y(N))
        ALLOCATE(z(N))
        ALLOCATE(forces0(3,N))
        ALLOCATE(forces(3,N))
        DO i_test = 1, tests + 1
            CALL generate_initial_values(N, m, x, y, z)

            ! compute forces directly
            ! CALL CPU_TIME(time_start)
            time_start = omp_get_wtime()
            ! forces0 = compute_forces_direct(m, x, y, z, 1.0_pf, eps)

            forces0 = compute_forces_direct(m, x, y, z, 1.0_pf, eps2, nt)

            ! !$OMP PARALLEL SHARED(forces0) PRIVATE(p) NUM_THREADS(nt)
            ! !$OMP DO
            ! DO p = 1, N
            !     forces0(:,p) = compute_forces_direct_over_p(m,x,y,z,1.0_pf,eps2,p)
            ! END DO
            ! !$OMP END DO
            ! !$OMP END PARALLEL

            ! CALL CPU_TIME(time_finish)
            time_finish = omp_get_wtime()
            total = time_finish - time_start
            IF (i_test > 1) WRITE(file, *) N, -1.0_pf, total, 0.0_pf

            ALLOCATE(tree)
            
            time_start = omp_get_wtime()
            CALL tree % init(m, x, y, z, multipole)
            time_finish = omp_get_wtime()
            time_generate_tree = time_finish - time_start

            DO i_theta = 1, SIZE(thetas)
                theta = thetas(i_theta)**2

                ! CALL CPU_TIME(time_start)
                time_start = omp_get_wtime()

                ! now test the tree
                IF (nt == 1) THEN
                    DO p = 1, N
                        IF (multipole == 1 .or. multipole == 4) THEN
                            forces(:,p) = tree % forces(p, theta, 1.0_pf, eps2)
                        ELSE
                            forces(:,p) = tree % evaluate_forces_over_p_octupole(p, theta, 1.0_pf, eps2)
                        ENDIF
                    END DO
                ELSE
                    !$OMP PARALLEL DO SHARED(forces) PRIVATE(p) NUM_THREADS(nt) &
                    !$OMP SCHEDULE(DYNAMIC)
                    DO p = 1, N
                        IF (multipole == 1 .or. multipole == 4) THEN
                            forces(:,p) = tree % forces(p, theta, 1.0_pf, eps2)
                        ELSE
                            forces(:,p) = tree % evaluate_forces_over_p_octupole(p, theta, 1.0_pf, eps2)
                        ENDIF
                    END DO
                    !$OMP END PARALLEL DO
                ENDIF

                ! CALL CPU_TIME(time_finish)
                time_finish = omp_get_wtime()

                total = time_finish - time_start + time_generate_tree

                ! erro = NORM2(forces - forces0)/NORM2(forces0)
                a_scale = 0.0_pf
                DO p = 1, N
                    a_scale = a_scale + m(p) / (x(p)**2 + y(p)**2 + z(p)**2)
                END DO
                erro = SQRT((NORM2(forces - forces0)/a_scale)**2 / N)

                IF (i_test > 1) WRITE(file, *) N, thetas(i_theta), total, erro
            END DO

            DEALLOCATE(tree)
        END DO
        DEALLOCATE(m, x, y, z, forces, forces0)
    END DO
END SUBROUTINE

FUNCTION compute_forces_direct_over_p (m, x, y, z, G, eps2, p) RESULT (forces)
    REAL(pf), INTENT(IN) :: m(:), x(:), y(:), z(:)
    REAL(pf), INTENT(IN) :: G, eps2
    INTEGER, INTENT(IN) :: p
    INTEGER :: b
    REAL(pf) :: dx, dy, dz, dist2, dist, f
    REAL(pf) :: forces(3)

    forces = 0.0_pf
    DO b = 1, SIZE(m)
        IF (b == p) CYCLE
        dx = x(b) - x(p)
        dy = y(b) - y(p)
        dz = z(b) - z(p)
        dist2 = dx*dx + dy*dy + dz*dz + eps2
        dist = SQRT(dist2)
        f = G * m(p) * m(b) / (dist * dist2)
        forces(1) = forces(1) + f * dx
        forces(2) = forces(2) + f * dy
        forces(3) = forces(3) + f * dz
    END DO
END FUNCTION

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

END PROGRAM
