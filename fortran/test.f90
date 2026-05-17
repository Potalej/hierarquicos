! # TEST
! 
! A simple test to generate the quadtree and save its information
! as a text file.
PROGRAM main

    USE tree_mod
    IMPLICIT NONE

    INTEGER, PARAMETER :: N = 30
    REAL(8) :: m(N), qs(N,2), x(N), y(N)
    INTEGER :: i, file = 45
    TYPE(TreeType) :: tree

    CALL RANDOM_SEED()

    m = 1.0d0

    CALL RANDOM_NUMBER(qs)
    qs = 2.0d0 * qs - 1.0d0

    ! separando colunas
    x = qs(:, 1)
    y = qs(:, 2)

    print *, "x =", x
    print *, "y =", y

    OPEN(file, file = "out/quadtree_test.txt", status="replace")

    WRITE(file, *) N
    WRITE(file, *) m
    WRITE(file, *) x
    WRITE(file, *) y

    CALL tree % init(m, x, y, file)
END PROGRAM