PROGRAM forces_time_and_error
    USE test_forces_mod
    IMPLICIT NONE

    INTEGER, PARAMETER :: pf  = SELECTED_REAL_KIND(15, 307)
    INTEGER :: file = 13

    CHARACTER(LEN=100) :: out_dir
    REAL(pf), ALLOCATABLE :: thetas(:)
    REAL(pf) :: eps2
    INTEGER :: nt, tests, ntidx
    INTEGER, ALLOCATABLE :: Nlist(:), nums_threads(:)
    LOGICAL :: eval_df
    
    CHARACTER(LEN=14) :: stamp

    CALL time_now(stamp)
    print *, "starting at ", stamp

    ALLOCATE(nums_threads(7))
    nums_threads = [1, 2, 4, 8, 12, 16, 20]

    tests = 6
    ALLOCATE(thetas(4))
    thetas = [0.25_pf, 0.5_pf, 0.75_pf, 1.0_pf]

    eps2 = 0.0_pf**2
    eval_df = .FALSE.

    ! CALL uniform_list(1000, 10000, 500, Nlist)
    CALL log_list(10000, 1000000, 18, Nlist)
    print *, "Nlist:", Nlist

    ! first we run the dehnen method because its only sequential
    out_dir = "out_brucutuvi/forces_only_time_dehnen/"
    CALL SYSTEM("mkdir -p "//TRIM(out_dir))

    CALL time_now(stamp)
    print *, '> DEHNEN METHOD ', stamp
    print *, '# monopole ', stamp
    OPEN(file, file = TRIM(out_dir)//"dehnen_test_multipoles_thetas_monopole.txt", status="replace")
        CALL test_forces_time(Nlist, thetas, eps2, tests, file, nt, 20, eval_df, 10)
    CLOSE(file)

    CALL time_now(stamp)
    print *, '# quadrupole ', stamp
    OPEN(file, file = TRIM(out_dir)//"dehnen_test_multipoles_thetas_quadrupole.txt", status="replace")
        CALL test_forces_time(Nlist, thetas, eps2, tests, file, nt, 21, eval_df, 10)
    CLOSE(file)

    CALL time_now(stamp)
    print *, '# octupole ', stamp
    OPEN(file, file = TRIM(out_dir)//"dehnen_test_multipoles_thetas_octupole.txt", status="replace")
        CALL test_forces_time(Nlist, thetas, eps2, tests, file, nt, 22, eval_df, 10)
    CLOSE(file)

    DO ntidx = 1, SIZE(nums_threads)
        nt = ntidx

        WRITE(out_dir, "(A,I0,A)") "out_brucutuvi/forces_only_time_p",nt,"/"
        CALL SYSTEM("mkdir -p "//TRIM(out_dir))

        ! now we test the direct forces
        CALL time_now(stamp)
        WRITE(*,"(A,I0,A,A)") '> DIRECT FORCES (', nt, ' threads) ', stamp
        OPEN(file, file = TRIM(out_dir)//"direct.txt", status="replace")
            CALL test_forces_time(Nlist, thetas, eps2, tests, file, nt, 0, .TRUE., nt)
        CLOSE(file)

        CALL time_now(stamp)
        WRITE(*,"(A,I0,A,A)") '> CLASSICAL BARNES-HUT METHOD (', nt, ' threads) ', stamp
            print *, '# monopole ', stamp
            OPEN(file, file = TRIM(out_dir)//"bh_test_multipoles_thetas_monopole.txt", status="replace")
                CALL test_forces_time(Nlist, thetas, eps2, tests, file, nt, 10, eval_df, 10)
            CLOSE(file)

            CALL time_now(stamp)
            print *, '# quadrupole ', stamp
            OPEN(file, file = TRIM(out_dir)//"bh_test_multipoles_thetas_quadrupole.txt", status="replace")
                CALL test_forces_time(Nlist, thetas, eps2, tests, file, nt, 11, eval_df, 10)
            CLOSE(file)

            CALL time_now(stamp)
            print *, '# octupole ', stamp
            OPEN(file, file = TRIM(out_dir)//"bh_test_multipoles_thetas_octupole.txt", status="replace")
                CALL test_forces_time(Nlist, thetas, eps2, tests, file, nt, 12, eval_df, 10)
            CLOSE(file)

        CALL time_now(stamp)
        WRITE(*,"(A,I0,A,A)") '> BARNES-HUT GROUPS METHOD (', nt, ' threads) ', stamp
            print *, '# monopole ', stamp
            OPEN(file, file = TRIM(out_dir)//"bhg_test_multipoles_thetas_monopole.txt", status="replace")
                CALL test_forces_time(Nlist, thetas, eps2, tests, file, nt, 30, eval_df, 10)
            CLOSE(file)

            CALL time_now(stamp)
            print *, '# quadrupole ', stamp
            OPEN(file, file = TRIM(out_dir)//"bhg_test_multipoles_thetas_quadrupole.txt", status="replace")
                CALL test_forces_time(Nlist, thetas, eps2, tests, file, nt, 31, eval_df, 10)
            CLOSE(file)

            CALL time_now(stamp)
            print *, '# octupole ', stamp
            OPEN(file, file = TRIM(out_dir)//"bhg_test_multipoles_thetas_octupole.txt", status="replace")
                CALL test_forces_time(Nlist, thetas, eps2, tests, file, nt, 32, eval_df, 10)
            CLOSE(file)
    END DO

CONTAINS

SUBROUTINE time_now (stamp)
    CHARACTER(LEN=14), INTENT(OUT) :: stamp
    INTEGER :: v(8)

    CALL DATE_AND_TIME(VALUES=v)
    WRITE(stamp, '(I2.2,"/",I2.2,"/",I2.2," ",I2.2,":",I2.2)') v(3), v(2), MOD(v(1), 100), v(5), v(6)
END SUBROUTINE

END PROGRAM