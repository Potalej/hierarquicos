PROGRAM forces_time_and_error
    USE octree_mod
    USE omp_lib
    USE initial_values_mod
    USE forces_mod
    IMPLICIT NONE

    INTEGER, PARAMETER :: pf  = SELECTED_REAL_KIND(15, 307)

    CALL compare_forces_methods()
CONTAINS

SUBROUTINE compare_forces_methods ()
    CHARACTER(LEN=100) :: out_dir
    REAL(pf) :: thetas(5), eps
    INTEGER  :: number_of_threads, number_of_tests
    INTEGER  :: Nmin, Nmax, Nstep
    INTEGER  :: file = 13
    
    NAMELIST /preset/ number_of_threads, number_of_tests, thetas, eps, Nmin, Nmax, Nstep

    out_dir = "out/direct_vs_bh_vs_dehnen/"
    CALL SYSTEM("mkdir -p "//TRIM(out_dir))

    number_of_threads = 1
    number_of_tests = 5
    eps = 0.0_pf
    thetas = (/0.2_pf, 0.4_pf, 0.6_pf, 0.8_pf, 1.0_pf/)

    Nmin = 5000
    Nmax = 50000
    Nstep = 5000

    ! preset file
    OPEN(file, file=TRIM(out_dir)//"preset.nml")
    WRITE(file, nml=preset)
    CLOSE(file)

    ! tests
    print *, '# bh monopole'
    OPEN(file, file = TRIM(out_dir)//"bh_monopole.txt", status="replace")
    CALL test_forces_time(Nmin, Nmax, Nstep, thetas, eps**2, number_of_tests, file, number_of_threads, 10)
    CLOSE(file)

    print *, '# bh quadrupole'
    OPEN(file, file = TRIM(out_dir)//"bh_quadrupole.txt", status="replace")
    CALL test_forces_time(Nmin, Nmax, Nstep, thetas, eps**2, number_of_tests, file, number_of_threads, 11)
    CLOSE(file)

    print *, '# dehnen monopole'
    OPEN(file, file = TRIM(out_dir)//"dehnen_monopole.txt", status="replace")
    CALL test_forces_time(Nmin, Nmax, Nstep, thetas, eps**2, number_of_tests, file, number_of_threads, 20)
    CLOSE(file)

    print *, '# dehnen quadrupole'
    OPEN(file, file = TRIM(out_dir)//"dehnen_quadrupole.txt", status="replace")
    CALL test_forces_time(Nmin, Nmax, Nstep, thetas, eps**2, number_of_tests, file, number_of_threads, 21)
    CLOSE(file)
END SUBROUTINE

!===============================================
! UTILITIES
!===============================================
SUBROUTINE test_forces_time (Nmin, Nmax, Nstep, thetas, eps2, tests, file, nt, method)
    INTEGER,  INTENT(IN) :: Nmin, Nmax, Nstep, tests, file, nt, method
    REAL(pf), INTENT(IN) :: thetas(:), eps2
    INTEGER :: N, i_test, i_theta, p
    REAL(pf), ALLOCATABLE :: m(:), qs(:,:), ps(:,:)       ! state vectors
    REAL(pf), ALLOCATABLE :: forces(:,:), forces_dir(:,:) ! forces

    REAL(pf) :: time_start, time_finish, total, time_generate_tree ! timers
    REAL(pf) :: theta, error, error2, time_quad, time_forces
    CLASS(OctreeType), ALLOCATABLE :: tree
    REAL(pf) :: a_scale

    DO N = Nmin, Nmax, Nstep
        PRINT *, 'N=', N
        ALLOCATE(m(N))
        ALLOCATE(qs(3,N), ps(3,N))
        ALLOCATE(forces(3,N), forces_dir(3,N))

        DO i_test = 1, tests + 1
            CALL generate_initial_values(N, m, qs, ps, &
                10.0_pf, & ! size_qs
                0.0_pf,  & ! size_ps
                .TRUE.,  & ! normalized
                .FALSE., & ! 2d
                "plummer", & ! mass density profile
                (/1.0_pf, 0.58_pf /) & ! profile parameters
            )

            ! compute forces directly
            time_start = omp_get_wtime()
                forces_dir = compute_forces_direct(m, qs(1,:), qs(2,:), qs(3,:), 1.0_pf, eps2, nt)
            time_finish = omp_get_wtime()
            total = time_finish - time_start
            WRITE(file, *) N, -1.0_pf, 0.0_pf, total, 0.0_pf, total, 0.0_pf

            ALLOCATE(tree)
            CALL tree % pre_init(m, method)

            time_start = omp_get_wtime()
                CALL tree % init(qs(1,:), qs(2,:), qs(3,:))
            time_finish = omp_get_wtime()
            time_generate_tree = time_finish - time_start

            DO i_theta = 1, SIZE(thetas)
                theta = thetas(i_theta)**2

                IF (method == 10 .OR. method == 11 .OR. method == 12) THEN
                    time_start = omp_get_wtime()
                        ! now test the tree
                        IF (nt == 1) THEN
                            DO p = 1, N
                                forces(:,p) = tree % forces(p, theta, 1.0_pf, eps2)
                            END DO
                        ELSE
                            !$OMP PARALLEL DO SHARED(forces) PRIVATE(p) NUM_THREADS(nt) &
                            !$OMP SCHEDULE(DYNAMIC)
                            DO p = 1, N
                                forces(:,p) = tree % forces(p, theta, 1.0_pf, eps2)
                            END DO
                            !$OMP END PARALLEL DO
                        ENDIF
                    time_finish = omp_get_wtime()
                    time_forces = time_finish - time_start
                    total = time_generate_tree + time_forces
                ELSE IF (method == 20 .OR. method == 21) THEN
                    tree % ns_force = 0.0_pf
                    IF (method == 21) THEN
                        tree % ns_hess = 0.0_pf
                        tree % ns_third = 0.0_pf
                    ENDIF
                    time_start = omp_get_wtime()
                        CALL tree % dehnen_eval(theta, eps2, 1.0_pf, forces)
                    time_finish = omp_get_wtime()
                    time_forces = time_finish - time_start
                    total = time_generate_tree + time_forces
                ELSE
                    STOP "method not identified!"
                ENDIF

                CALL evaluate_error_on_accelerations(m, qs, forces, forces_dir, error, error2)

                WRITE(file, *) N, thetas(i_theta), time_generate_tree, total, error, time_forces, error2
            END DO

            DEALLOCATE(tree)
        END DO
        DEALLOCATE(m, qs, ps, forces, forces_dir)
    END DO
END SUBROUTINE

SUBROUTINE evaluate_error_on_accelerations (ms, qs, forces, forces_dir, error, error2)
    REAL(pf), INTENT(IN)  :: ms(:), qs(:,:)
    REAL(pf), INTENT(OUT) :: error, error2
    REAL(pf), INTENT(IN)  :: forces(:,:), forces_dir(:,:)
    REAL(pf) :: accel(3, SIZE(ms)), accel_dir(3, SIZE(ms))
    INTEGER  :: p, N
    REAL(pf) :: a_scale

    N = size(ms)
    
    a_scale = 0.0_pf
    DO p = 1, N
        ! transforming to accelerations
        accel(:,p) = forces(:,p) / ms(p)
        accel_dir(:,p) = forces_dir(:,p) / ms(p)

        ! scale
        a_scale = a_scale + ms(p) / DOT_PRODUCT(qs(:,p), qs(:,p))
    END DO

    error = SQRT((NORM2(accel - accel_dir)/a_scale)**2 / N)

    a_scale = 0.0_pf
    DO p = 1, N
        ! transforming to accelerations
        accel(:,p) = forces(:,p) / ms(p)
        accel_dir(:,p) = forces_dir(:,p) / ms(p)
    END DO
    error2 = NORM2(accel - accel_dir)/NORM2(accel_dir)
END SUBROUTINE

END PROGRAM