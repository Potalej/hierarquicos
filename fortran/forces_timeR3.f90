! # TREE TIME
! 
! Program to test the time to construct the tree
! with N in [10, 1e4] (10 by 10) with 50 tests
! by value of N.
! It achieves the O(N log N), as expected :)
PROGRAM forces_time
    USE octree_mod
    USE morton_octree_mod
    USE omp_lib
    IMPLICIT NONE

    INTEGER :: file
    REAL(8) :: thetas(4)

    file = 45

    thetas(1) = 0.0d0
    thetas(2) = 0.2d0
    thetas(3) = 0.25d0
    thetas(4) = 0.5d0

    !OPEN(file, file = "out/forces_time_morton_parallel_20C_N10K.txt", status="replace")
    !CALL test_forces_time_morton(500, 10000, 500, thetas, 0.1d0, 20, file)
    !CLOSE(file)

    OPEN(file, file = "out/forces_time_without_morton_parallel_20C_N10K.txt", status="replace")
    CALL test_forces_time(100, 20000, 100, thetas, 0.1d0, 10, file)
    CLOSE(file)
CONTAINS

SUBROUTINE generate_initial_values (N, m, x, y, z)
    INTEGER, INTENT(IN) :: N
    REAL(8), INTENT(INOUT) :: m(N), x(N), y(N), z(N)
    REAL(8) :: qs(N,3)

    CALL RANDOM_SEED()
    m = 1.0d0

    CALL RANDOM_NUMBER(qs)
    qs = 10.0d0 * (2.0d0 * qs - 1.0d0)

    x = qs(:, 1)
    y = qs(:, 2)
    z = qs(:, 3)
END SUBROUTINE

SUBROUTINE test_forces_time (Nmin, Nmax, Nstep, thetas, eps, tests, file)
    INTEGER, INTENT(IN) :: Nmin, Nmax, Nstep, tests, file
    REAL(8), INTENT(IN) :: thetas(:), eps
    INTEGER :: N, i_test, i_theta, p, timer
    REAL(8), ALLOCATABLE :: m(:), x(:), y(:), z(:), forces(:,:), forces0(:,:)
    REAL(8) :: time_start, time_finish, total
    REAL(8) :: theta, erro, time_quad
    CLASS(OctreeType), ALLOCATABLE :: tree

    DO N = Nmin, Nmax, Nstep
        PRINT *, 'N=', N
        ALLOCATE(m(N))
        ALLOCATE(x(N))
        ALLOCATE(y(N))
        ALLOCATE(z(N))
        ALLOCATE(forces0(N,3))
        ALLOCATE(forces(N,3))
        DO i_test = 1, tests + 1
            CALL generate_initial_values(N, m, x, y, z)

            ! compute forces directly
            ! CALL CPU_TIME(time_start)
            time_start = omp_get_wtime()
            ! forces0 = compute_forces_direct(m, x, y, z, 1.0d0, eps)

            !$OMP PARALLEL SHARED(forces0) PRIVATE(p) NUM_THREADS(20)
            !$OMP DO
            DO p = 1, N
                forces0(p,:) = compute_forces_direct_over_p(m,x,y,z,1.0d0,eps,p)
            END DO
            !$OMP END DO
            !$OMP END PARALLEL

            ! CALL CPU_TIME(time_finish)
            time_finish = omp_get_wtime()
            total = time_finish - time_start
            IF (i_test > 1) WRITE(file, *) N, -1.0d0, total, 0.0d0

            ALLOCATE(tree)
            
            time_start = omp_get_wtime()
            CALL tree % init(m, x, y, z)
            time_finish = omp_get_wtime()

            ! print *, ''
            ! print *, 'init:', time_finish - time_start

            time_start = omp_get_wtime()
            CALL tree % evaluate_quad()
            time_finish = omp_get_wtime()
            ! print *, 'quad:', time_finish - time_start

            DO i_theta = 1, SIZE(thetas)
                theta = thetas(i_theta)**2

                ! CALL CPU_TIME(time_start)
                time_start = omp_get_wtime()

                ! now test the tree
                !$OMP PARALLEL SHARED(forces) PRIVATE(p) NUM_THREADS(20)
                !$OMP DO
                DO p = 1, N
                    forces(p,:) = tree % forces(p, theta, 1.0d0, eps)
                END DO
                !$OMP END DO
                !$OMP END PARALLEL

                ! CALL CPU_TIME(time_finish)
                time_finish = omp_get_wtime()

                total = time_finish - time_start

                erro = NORM2(forces - forces0)/NORM2(forces0)

                IF (i_test > 1) WRITE(file, *) N, thetas(i_theta), total, erro
            END DO

            DEALLOCATE(tree)
        END DO
        DEALLOCATE(m, x, y, z, forces, forces0)
    END DO
END SUBROUTINE

SUBROUTINE test_forces_time_morton (Nmin, Nmax, Nstep, thetas, eps, tests, file)
    INTEGER, INTENT(IN) :: Nmin, Nmax, Nstep, tests, file
    REAL(8), INTENT(IN) :: thetas(:), eps
    INTEGER :: N, i_test, i_theta, p, timer
    REAL(8), ALLOCATABLE :: m(:), x(:), y(:), z(:), forces(:,:), forces0(:,:), qs(:,:)
    REAL(8) :: time_start, time_finish, total
    REAL(8) :: theta, erro, time_quad
    INTEGER :: L = 0

    TYPE(morton_tree_type), ALLOCATABLE :: morton_tree

    DO N = Nmin, Nmax, Nstep
        PRINT *, 'N=', N
        ALLOCATE(qs(N,3))
        ALLOCATE(m(N))
        ALLOCATE(x(N))
        ALLOCATE(y(N))
        ALLOCATE(z(N))
        ALLOCATE(forces0(N,3))
        ALLOCATE(forces(N,3))

        DO i_test = 1, tests + 1
            CALL generate_initial_values(N, m, x, y, z)
            qs(:,1) = x
            qs(:,2) = y
            qs(:,3) = z

            ! compute forces directly
            ! CALL CPU_TIME(time_start)
            time_start = omp_get_wtime()

            !$OMP PARALLEL SHARED(forces) PRIVATE(p) NUM_THREADS(20)
            !$OMP DO SCHEDULE(STATIC)
            DO p = 1, N
                forces0(p,:) = compute_forces_direct_over_p(m,x,y,z,1.0d0,eps,p)
            END DO
            !$OMP END DO
            !$OMP END PARALLEL

            ! CALL CPU_TIME(time_finish)
            time_finish = omp_get_wtime()
            total = time_finish - time_start
            IF (i_test > 1) WRITE(file, *) N, -1.0d0, total, 0.0d0

            ALLOCATE(morton_tree)
            CALL morton_tree % init(N, m, qs, L, 1.2d0, .true.)

            DO i_theta = 1, SIZE(thetas)
                theta = thetas(i_theta)**2

                ! CALL CPU_TIME(time_start)
                time_start = omp_get_wtime()

                forces = morton_tree % forces_par(1.0d0, eps**2, theta, .true., 20)

                ! CALL CPU_TIME(time_finish)
                time_finish = omp_get_wtime()

                total = time_finish - time_start

                erro = NORM2(forces - forces0)/NORM2(forces0)

                IF (i_test > 1) WRITE(file, *) N, thetas(i_theta), total, erro
            END DO

            DEALLOCATE(morton_tree)
        END DO
        DEALLOCATE(m, qs, x, y, z, forces, forces0)
    END DO
END SUBROUTINE

FUNCTION compute_forces_direct_over_p (m, x, y, z, G, eps, p) RESULT (forces)
    REAL(8), INTENT(IN) :: m(:), x(:), y(:), z(:)
    REAL(8), INTENT(IN) :: G, eps
    INTEGER, INTENT(IN) :: p
    INTEGER :: b
    REAL(8) :: dx, dy, dz, dist2, dist, f
    REAL(8) :: forces(3)

    forces = 0.0d0
    DO b = 1, SIZE(m)
        IF (b == p) CYCLE
        dx = x(b) - x(p)
        dy = y(b) - y(p)
        dz = z(b) - z(p)
        dist2 = dx*dx + dy*dy + dz*dz + eps*eps
        dist = SQRT(dist2)
        f = G * m(p) * m(b) / (dist * dist2)
        forces(1) = forces(1) + f * dx
        forces(2) = forces(2) + f * dy
        forces(3) = forces(3) + f * dz
    END DO
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

END PROGRAM
