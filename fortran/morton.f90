! morton sort. im not sure how this really works
MODULE morton_mod
IMPLICIT NONE
INTEGER, PARAMETER :: i8 = selected_int_kind(18)

CONTAINS

PURE FUNCTION split_by_3(a) RESULT(x)
    INTEGER(i8), INTENT(IN) :: a
    INTEGER(i8) :: x

    x = a

    x = IAND(x, INT(Z'1FFFFF',i8))

    x = IAND(IOR(x, ISHFT(x,32)), INT(Z'1F00000000FFFF',i8))
    x = IAND(IOR(x, ISHFT(x,16)), INT(Z'1F0000FF0000FF',i8))
    x = IAND(IOR(x, ISHFT(x, 8)), INT(Z'100F00F00F00F00F',i8))
    x = IAND(IOR(x, ISHFT(x, 4)), INT(Z'10C30C30C30C30C3',i8))
    x = IAND(IOR(x, ISHFT(x, 2)), INT(Z'1249249249249249',i8))
END FUNCTION

PURE FUNCTION morton3D(ix, iy, iz) RESULT(code)
    INTEGER(i8), INTENT(IN) :: ix, iy, iz
    INTEGER(i8) :: code

    code = split_by_3(ix)                      &
        + ISHFT(split_by_3(iy), 1)            &
        + ISHFT(split_by_3(iz), 2)
END FUNCTION

SUBROUTINE morton_sort_3d(x, y, z, m)

    REAL(8), INTENT(INOUT) :: x(:), y(:), z(:), m(:)
    INTEGER :: N
    INTEGER :: p
    INTEGER(i8), ALLOCATABLE :: key(:)
    INTEGER, ALLOCATABLE :: idx(:)
    REAL(8), ALLOCATABLE :: xt(:), yt(:), zt(:), mt(:)
    REAL(8) :: xmin,xmax
    REAL(8) :: ymin,ymax
    REAL(8) :: zmin,zmax
    REAL(8) :: scale
    INTEGER(i8) :: ix,iy,iz

    N = SIZE(x)

    ALLOCATE(key(N))
    ALLOCATE(idx(N))

    !-----------------------------------------
    ! normalize coordinates
    !-----------------------------------------

    xmin = MINVAL(x)
    xmax = MAXVAL(x)

    ymin = MINVAL(y)
    ymax = MAXVAL(y)

    zmin = MINVAL(z)
    zmax = MAXVAL(z)

    scale = REAL(2_i8**21 - 1,8)

    !-----------------------------------------
    ! compute morton keys
    !-----------------------------------------

    DO p=1,N

        ix = INT(scale * (x(p)-xmin)/(xmax-xmin), i8)
        iy = INT(scale * (y(p)-ymin)/(ymax-ymin), i8)
        iz = INT(scale * (z(p)-zmin)/(zmax-zmin), i8)

        key(p) = morton3D(ix,iy,iz)

        idx(p) = p

    ENDDO

    !-----------------------------------------
    ! sort indices by key
    !-----------------------------------------

    CALL quicksort(key, idx, 1, N)

    !-----------------------------------------
    ! permute arrays
    !-----------------------------------------

    ALLOCATE(xt(N), yt(N), zt(N), mt(N))

    xt = x(idx)
    yt = y(idx)
    zt = z(idx)
    mt = m(idx)

    x = xt
    y = yt
    z = zt
    m = mt

    DEALLOCATE(key, idx)
    DEALLOCATE(xt, yt, zt, mt)
END SUBROUTINE

RECURSIVE SUBROUTINE quicksort(key, idx, left, right)

    INTEGER(i8), INTENT(INOUT) :: key(:)
    INTEGER, INTENT(INOUT) :: idx(:)

    INTEGER, INTENT(IN) :: left, right

    INTEGER :: i, j
    INTEGER(i8) :: pivot, tempk
    INTEGER :: tempi

    IF (left >= right) RETURN

    pivot = key((left+right)/2)

    i = left
    j = right

    DO

        DO WHILE (key(i) < pivot)
            i = i + 1
        ENDDO

        DO WHILE (key(j) > pivot)
            j = j - 1
        ENDDO

        IF (i <= j) THEN

            tempk = key(i)
            key(i) = key(j)
            key(j) = tempk

            tempi = idx(i)
            idx(i) = idx(j)
            idx(j) = tempi

            i = i + 1
            j = j - 1

        ENDIF

        IF (i > j) EXIT

    ENDDO

    IF (left < j) CALL quicksort(key, idx, left, j)
    IF (i < right) CALL quicksort(key, idx, i, right)

END SUBROUTINE

END MODULE