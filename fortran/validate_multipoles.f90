! # VALIDATING THE MULTIPOLES IMPLEMENTATIONS
PROGRAM forces_time
    USE octree_mod
    USE omp_lib
    USE forces_mod
    IMPLICIT NONE

    INTEGER, PARAMETER :: pf  = SELECTED_REAL_KIND(15, 307)

    INTEGER :: file, nt
    REAL(pf) :: thetas(9), eps2

    INTEGER :: N
    LOGICAL :: quad

    file = 45
    nt = 10 ! number of threads
    ! eps2 = (0.1_pf)**2
    eps2 = 0.0_pf

    thetas(1) = 0.2_pf
    thetas(2) = 0.3_pf
    thetas(3) = 0.4_pf
    thetas(4) = 0.5_pf
    thetas(5) = 0.6_pf
    thetas(6) = 0.7_pf
    thetas(7) = 0.8_pf
    thetas(8) = 0.9_pf
    thetas(9) = 1.0_pf

    N = 10000

    CALL validate_multipoles_construction(N, thetas, eps2, 1)
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

SUBROUTINE validate_multipoles_construction (N, thetas, eps2, tests)
    INTEGER,  INTENT(IN) :: N, tests
    REAL(pf), INTENT(IN) :: eps2, thetas(:)
    INTEGER :: i_test, i_theta, timer
    REAL(pf) :: m(N), x(N), y(N), z(N)
    REAL(pf) :: time_start, time_finish, total, time_generate_tree
    REAL(pf) :: erro, time_quad
    CLASS(OctreeType), ALLOCATABLE :: tree
    REAL(pf) :: a_scale

    REAL(pf) :: forces_direct(3,N)
    REAL(pf) :: forces_0(3,N)
    REAL(pf) :: forces_1(3,N)
    REAL(pf) :: forces_4(3,N)
    REAL(pf) :: forces_8(3,N)
    INTEGER :: p

    PRINT *, 'N=', N

    CALL generate_initial_values(N, m, x, y, z)

    x(2) = x(1) + 1e-4
    y(2) = y(1) + 1e-4
    z(2) = z(1) + 1e-4

    ! direct forces
    forces_direct = compute_forces_direct(m, x, y, z, 1.0_pf, eps2, 1)

    ALLOCATE(tree)

    i_theta = 1

    ! tree monopole
    print *, "[MONOPOLE]"
    time_start = omp_get_wtime()
    CALL tree % init(m, x, y, z, 1)
    time_finish = omp_get_wtime()
    print *, 'time (construct)', time_finish - time_start

    time_start = omp_get_wtime()
    DO p = 1, N
        forces_1(:,p) = tree % forces(p, thetas(i_theta)**2, 1.0_pf, eps2)
    END DO
    time_finish = omp_get_wtime()
    print *, 'time (evaluate)', time_finish - time_start
    print *, "erros:", NORM2(forces_direct - forces_1)

    DEALLOCATE(tree)
    print *, ''
    ALLOCATE(tree)

    ! tree quadrupole
    print *, "[QUADRUPOLE]"
    time_start = omp_get_wtime()
    CALL tree % init(m, x, y, z, 4)
    time_finish = omp_get_wtime()
    print *, 'time (construct)', time_finish - time_start

    time_start = omp_get_wtime()
    DO p = 1, N
        forces_4(:,p) = tree % forces(p, thetas(i_theta)**2, 1.0_pf, eps2)
    END DO
    time_finish = omp_get_wtime()
    print *, 'time (evaluate)', time_finish - time_start
    print *, "erros:", NORM2(forces_direct - forces_4)

    DEALLOCATE(tree)
    print *, ''
    ALLOCATE(tree)

    ! tree octupole
    print *, "[OCTUPOLE]"
    time_start = omp_get_wtime()
    CALL tree % init(m, x, y, z, 8)
    time_finish = omp_get_wtime()
    print *, 'time (construct)', time_finish - time_start

    time_start = omp_get_wtime()
    DO p = 1, N
        forces_8(:,p) = tree % forces(p, thetas(i_theta)**2, 1.0_pf, eps2)
    END DO
    time_finish = omp_get_wtime()
    print *, 'time (evaluate)', time_finish - time_start
    print *, "erros:", NORM2(forces_direct - forces_8)

END SUBROUTINE


! SUBROUTINE validate_multipoles_construction (N, thetas, eps2, tests)
!     INTEGER,  INTENT(IN) :: N, tests
!     REAL(pf), INTENT(IN) :: eps2, thetas(:)
!     INTEGER :: i_test, i_theta, timer
!     REAL(pf) :: m(N), x(N), y(N), z(N), forces(3,N)
!     REAL(pf) :: time_start, time_finish, total, time_generate_tree
!     REAL(pf) :: erro, time_quad
!     CLASS(OctreeType), ALLOCATABLE :: tree_orig, tree_new
!     REAL(pf) :: a_scale

!     PRINT *, 'N=', N

!     CALL generate_initial_values(N, m, x, y, z)

!     ALLOCATE(tree_orig)
!     time_start = omp_get_wtime()
!     CALL tree_orig % init(m, x, y, z, 4, 1)
!     time_finish = omp_get_wtime()
!     print *, 'tempo original (quadrupolo):', time_finish - time_start

!     ALLOCATE(tree_new)
!     time_start = omp_get_wtime()
!     CALL tree_new % init(m, x, y, z, 4, 4)
!     time_finish = omp_get_wtime()
!     print *, 'tempo novo (quadrupolo):', time_finish - time_start

!     DEALLOCATE(tree_orig)
!     DEALLOCATE(tree_new)

!     ALLOCATE(tree_orig)
!     time_start = omp_get_wtime()
!     CALL tree_orig % init(m, x, y, z, 8, 1)
!     time_finish = omp_get_wtime()
!     print *, 'tempo original (octupolo):', time_finish - time_start

!     ALLOCATE(tree_new)
!     time_start = omp_get_wtime()
!     CALL tree_new % init(m, x, y, z, 8, 2)
!     time_finish = omp_get_wtime()
!     print *, 'tempo novo (octupolo):', time_finish - time_start
    
!     print *, 'norm2 quad:', NORM2(tree_orig % ns_quad(:,:) - tree_new % ns_quad(:,:))
!     print *, 'norm2 octu:', NORM2(tree_orig % ns_oct(:,:) - tree_new % ns_oct(:,:))

! END SUBROUTINE

END PROGRAM
