PROGRAM forces_time_and_error
    USE octree_mod
    USE omp_lib
    USE initial_values_mod
    USE forces_mod
    IMPLICIT NONE

    INTEGER, PARAMETER :: pf  = SELECTED_REAL_KIND(15, 307)

    ! Test: Reproducing the results of McMillan and Aarseth.
    CALL test_mcmillan_aarseth_results()

    ! Test: N between 1e3 and 5e4
    CALL test_multipoles_bigger_N()

    ! Test: Sequential
    CALL test_multipoles_sequential()
CONTAINS

SUBROUTINE test_mcmillan_aarseth_results ()
! Here I'll reproduce the results of Aarseths and McMillan in the O(NlogN) for
! collisional stellars systems paper.
! For this I'll test for N=1024 and N=4096 and some angles values between 0.1 
! and 1.0
    CHARACTER(LEN=100) :: out_dir
    REAL(pf) :: thetas(9), eps
    INTEGER  :: number_of_threads, number_of_tests
    INTEGER  :: Nmin, Nmax, Nstep
    INTEGER  :: file = 13
    
    NAMELIST /preset/ number_of_threads, number_of_tests, thetas, eps, Nmin, Nmax, Nstep

    out_dir = "out/aarseth/"
    CALL SYSTEM("mkdir -p "//TRIM(out_dir))

    number_of_threads = 20
    number_of_tests = 20
    eps = 0.0_pf
    thetas = (/0.2_pf, 0.3_pf, 0.4_pf, 0.5_pf, 0.6_pf, 0.7_pf, 0.8_pf, 0.9_pf, 1.0_pf/)

    Nmin = 1024
    Nmax = 4096
    Nstep = 4096 - 1024

    ! preset file
    OPEN(file, file=TRIM(out_dir)//"preset.nml")
    WRITE(file, nml=preset)
    CLOSE(file)

    ! tests
    print *, '# monopole'
    OPEN(file, file = TRIM(out_dir)//"test_monopole.txt", status="replace")
    CALL test_forces_time(Nmin, Nmax, Nstep, thetas, eps**2, number_of_tests, file, number_of_threads, 1)
    CLOSE(file)

    print *, '# quadrupole'
    OPEN(file, file = TRIM(out_dir)//"test_quadrupole.txt", status="replace")
    CALL test_forces_time(Nmin, Nmax, Nstep, thetas, eps**2, number_of_tests, file, number_of_threads, 4)
    CLOSE(file)

    print *, '# octupole'
    OPEN(file, file = TRIM(out_dir)//"test_octupole.txt", status="replace")
    CALL test_forces_time(Nmin, Nmax, Nstep, thetas, eps**2, number_of_tests, file, number_of_threads, 8)
    CLOSE(file)
END SUBROUTINE

SUBROUTINE test_multipoles_bigger_N ()
! Here I'll test the octree for a bigger N, going from 1e3 to 5e4 (parallel).
    CHARACTER(LEN=100) :: out_dir
    REAL(pf) :: thetas(4), eps
    INTEGER  :: number_of_threads, number_of_tests
    INTEGER  :: Nmin, Nmax, Nstep
    INTEGER  :: file = 13
    
    NAMELIST /preset/ number_of_threads, number_of_tests, thetas, eps, Nmin, Nmax, Nstep

    out_dir = "out/bigger_N/"
    CALL SYSTEM("mkdir -p "//TRIM(out_dir))

    number_of_threads = 25
    number_of_tests = 20
    eps = 0.0_pf
    thetas = (/0.25_pf, 0.5_pf, 0.75_pf, 1.0_pf/)

    Nmin = 1000
    Nmax = 50000
    Nstep = 1000

    ! preset file
    OPEN(file, file=TRIM(out_dir)//"preset.nml")
    WRITE(file, nml=preset)
    CLOSE(file)

    ! tests
    print *, '# monopole'
    OPEN(file, file = TRIM(out_dir)//"test_monopole.txt", status="replace")
    CALL test_forces_time(Nmin, Nmax, Nstep, thetas, eps**2, number_of_tests, file, number_of_threads, 1)
    CLOSE(file)

    print *, '# quadrupole'
    OPEN(file, file = TRIM(out_dir)//"test_quadrupole.txt", status="replace")
    CALL test_forces_time(Nmin, Nmax, Nstep, thetas, eps**2, number_of_tests, file, number_of_threads, 4)
    CLOSE(file)

    print *, '# octupole'
    OPEN(file, file = TRIM(out_dir)//"test_octupole.txt", status="replace")
    CALL test_forces_time(Nmin, Nmax, Nstep, thetas, eps**2, number_of_tests, file, number_of_threads, 8)
    CLOSE(file)
END SUBROUTINE

SUBROUTINE test_multipoles_sequential ()
! Here I'll test the octree like the McMillan and Aarseth's test, but sequentially.
    CHARACTER(LEN=100) :: out_dir
    REAL(pf) :: thetas(5), eps
    INTEGER  :: number_of_threads, number_of_tests
    INTEGER  :: Nmin, Nmax, Nstep
    INTEGER  :: file = 13
    
    NAMELIST /preset/ number_of_threads, number_of_tests, thetas, eps, Nmin, Nmax, Nstep

    out_dir = "out/sequential/"
    CALL SYSTEM("mkdir -p "//TRIM(out_dir))

    number_of_threads = 1
    number_of_tests = 20
    eps = 0.0_pf
    thetas = (/0.2_pf, 0.4_pf, 0.6_pf, 0.8_pf, 1.0_pf/)

    Nmin = 500
    Nmax = 10000
    Nstep = 500

    ! preset file
    OPEN(file, file=TRIM(out_dir)//"preset.nml")
    WRITE(file, nml=preset)
    CLOSE(file)

    ! tests
    print *, '# monopole'
    OPEN(file, file = TRIM(out_dir)//"test_monopole.txt", status="replace")
    CALL test_forces_time(Nmin, Nmax, Nstep, thetas, eps**2, number_of_tests, file, number_of_threads, 1)
    CLOSE(file)

    print *, '# quadrupole'
    OPEN(file, file = TRIM(out_dir)//"test_quadrupole.txt", status="replace")
    CALL test_forces_time(Nmin, Nmax, Nstep, thetas, eps**2, number_of_tests, file, number_of_threads, 4)
    CLOSE(file)

    print *, '# octupole'
    OPEN(file, file = TRIM(out_dir)//"test_octupole.txt", status="replace")
    CALL test_forces_time(Nmin, Nmax, Nstep, thetas, eps**2, number_of_tests, file, number_of_threads, 8)
    CLOSE(file)
END SUBROUTINE

!===============================================
! UTILITIES
!===============================================
SUBROUTINE test_forces_time (Nmin, Nmax, Nstep, thetas, eps2, tests, file, nt, multipole)
    INTEGER,  INTENT(IN) :: Nmin, Nmax, Nstep, tests, file, nt
    REAL(pf), INTENT(IN) :: thetas(:), eps2
    INTEGER,  INTENT(IN) :: multipole
    INTEGER :: N, i_test, i_theta, p
    REAL(pf), ALLOCATABLE :: m(:), qs(:,:), ps(:,:)       ! state vectors
    REAL(pf), ALLOCATABLE :: forces(:,:), forces_dir(:,:) ! forces

    REAL(pf) :: time_start, time_finish, total, time_generate_tree ! timers
    REAL(pf) :: theta, error, time_quad
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
            WRITE(file, *) N, -1.0_pf, 0.0_pf, total, 0.0_pf

            ALLOCATE(tree)

            time_start = omp_get_wtime()
                CALL tree % init(m, qs(1,:), qs(2,:), qs(3,:), multipole)
            time_finish = omp_get_wtime()
            time_generate_tree = time_finish - time_start

            DO i_theta = 1, SIZE(thetas)
                theta = thetas(i_theta)**2

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
                total = time_finish - time_start + time_generate_tree

                CALL evaluate_error_on_accelerations(m, qs, forces, forces_dir, error)

                WRITE(file, *) N, thetas(i_theta), time_generate_tree, total, error
            END DO

            DEALLOCATE(tree)
        END DO
        DEALLOCATE(m, qs, ps, forces, forces_dir)
    END DO
END SUBROUTINE

SUBROUTINE evaluate_error_on_accelerations (ms, qs, forces, forces_dir, error)
    REAL(pf), INTENT(IN)  :: ms(:), qs(:,:)
    REAL(pf), INTENT(OUT) :: error
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
END SUBROUTINE

END PROGRAM