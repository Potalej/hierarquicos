! # TREE TIME
! 
! Program to test the time to construct the tree
! with N in [10, 1e4] (10 by 10) with 50 tests
! by value of N.
! It achieves the O(N log N), as expected :)
PROGRAM tree_time
    USE tree_mod
    IMPLICIT NONE

    INTEGER :: file

    file = 45

    OPEN(file, file = "out/tree_time.txt", status="replace")
    CALL test_tree_time(10, 10000, 10, 50, file)
CONTAINS

SUBROUTINE generate_initial_values (N, m, x, y)
    INTEGER, INTENT(IN) :: N
    REAL(8), INTENT(INOUT) :: m(N), x(N), y(N)
    REAL(8) :: qs(N,2)

    CALL RANDOM_SEED()
    m = 1.0d0

    CALL RANDOM_NUMBER(qs)
    qs = 2.0d0 * qs - 1.0d0

    x = qs(:, 1)
    y = qs(:, 2)
END SUBROUTINE

SUBROUTINE test_tree_time (Nmin, Nmax, Nstep, tests, file)
    INTEGER, INTENT(IN) :: Nmin, Nmax, Nstep, tests, file
    INTEGER :: N, i
    REAL(8), ALLOCATABLE :: m(:), x(:), y(:)
    REAL :: time_start, time_finish, total
    CLASS(TreeType), ALLOCATABLE :: tree

    DO N = Nmin, Nmax, Nstep
        DO i = 1, tests
            ALLOCATE(m(N))
            ALLOCATE(x(N))
            ALLOCATE(y(N))
            
            CALL generate_initial_values(N, m, x, y)

            CALL CPU_TIME(time_start)

            ! now test the tree
            ALLOCATE(tree)
            CALL tree % init(m, x, y)

            CALL CPU_TIME(time_finish)

            total = time_finish - time_start

            WRITE(file, *) N, total

            DEALLOCATE(tree)
            DEALLOCATE(m)
            DEALLOCATE(x)
            DEALLOCATE(y)
        END DO
    END DO
END SUBROUTINE

END PROGRAM