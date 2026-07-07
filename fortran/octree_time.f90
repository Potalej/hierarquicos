! # TREE TIME
! 
! Program to test the time to construct the tree
! with N in [10, 1e4] (10 by 10) with 50 tests
! by value of N.
! It achieves the O(N log N), as expected :)
PROGRAM octree_time
    USE octree_mod
    USE morton_octree_mod
    USE omp_lib
    IMPLICIT NONE

    INTEGER :: file

    file = 45

    OPEN(file, file = "out/octree_time_morton.txt", status="replace")
    ! CALL test_tree_time(1000, 20000, 100, 1, file)
    CALL test_tree_time_morton(1000, 20000, 100, 1, file)
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

SUBROUTINE test_tree_time (Nmin, Nmax, Nstep, tests, file)
    INTEGER, INTENT(IN) :: Nmin, Nmax, Nstep, tests, file
    INTEGER :: N, i
    REAL(8), ALLOCATABLE :: m(:), x(:), y(:), z(:)
    REAL :: time_start, time_finish, total
    CLASS(OctreeType), ALLOCATABLE :: tree

    DO N = Nmin, Nmax, Nstep
        DO i = 1, tests
            ALLOCATE(m(N))
            ALLOCATE(x(N))
            ALLOCATE(y(N))
            ALLOCATE(z(N))
            
            CALL generate_initial_values(N, m, x, y, z)

            CALL CPU_TIME(time_start)

            ! now test the tree
            ALLOCATE(tree)
            CALL tree % init(m, x, y, z)
            CALL tree % evaluate_quad()

            CALL CPU_TIME(time_finish)

            total = time_finish - time_start

            WRITE(file, *) N, total
            print *, N, total

            DEALLOCATE(tree)
            DEALLOCATE(m)
            DEALLOCATE(x, y, z)
        END DO
    END DO
END SUBROUTINE

SUBROUTINE test_tree_time_morton (Nmin, Nmax, Nstep, tests, file)
    INTEGER, INTENT(IN) :: Nmin, Nmax, Nstep, tests, file
    INTEGER :: N, i
    REAL(8), ALLOCATABLE :: m(:), x(:), y(:), z(:)
    REAL :: time_start, time_finish, total
    TYPE(node_type), ALLOCATABLE :: nodes(:)
    INTEGER :: nnodes
    REAL(8), ALLOCATABLE :: qs(:,:)
    INTEGER :: L

    TYPE(morton_tree_type), ALLOCATABLE :: morton_tree

    L = 0

    DO N = Nmin, Nmax, Nstep
        DO i = 1, tests
            ALLOCATE(m(N))
            ALLOCATE(x(N))
            ALLOCATE(y(N))
            ALLOCATE(z(N))
            ALLOCATE(qs(N,3))
            
            CALL generate_initial_values(N, m, x, y, z)
            qs(:,1) = x
            qs(:,2) = y
            qs(:,3) = z

            CALL CPU_TIME(time_start)

            ! now test the tree
            ALLOCATE(morton_tree)
            CALL morton_tree % init(N, m, qs, L, 1.2d0, .true.)

            CALL CPU_TIME(time_finish)

            total = time_finish - time_start

            WRITE(file, *) N, total
            PRINT *, N, total

            DEALLOCATE(m)
            DEALLOCATE(x, y, z, qs)
            DEALLOCATE(morton_tree)
        END DO
    END DO
END SUBROUTINE

END PROGRAM