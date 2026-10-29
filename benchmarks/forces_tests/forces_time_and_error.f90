PROGRAM forces_time_and_error
    USE test_forces_mod
    IMPLICIT NONE

    INTEGER, PARAMETER :: pf  = SELECTED_REAL_KIND(15, 307)
    INTEGER :: file = 13

    CHARACTER(LEN=100) :: out_dir
    REAL(pf), ALLOCATABLE :: thetas(:)
    REAL(pf) :: eps2
    INTEGER :: nt, tests
    INTEGER, ALLOCATABLE :: Nlist(:)
    LOGICAL :: eval_df

    out_dir = "out/forces_time_and_error/"
    CALL SYSTEM("mkdir -p "//TRIM(out_dir))

    nt = 20
    eps2 = 0.0_pf**2
    eval_df = .TRUE.

    tests = 5
    ALLOCATE(thetas(4))
    thetas = [0.25_pf, 0.5_pf, 0.75_pf, 1.0_pf]

    ! CALL uniform_list(1000, 10000, 500, Nlist)
    CALL log_list(1000, 100000, 14, Nlist)
    print *, "Nlist:", Nlist

    print *, '> CLASSICAL BARNES-HUT METHOD'
        print *, '# monopole'
        OPEN(file, file = TRIM(out_dir)//"bh_test_multipoles_thetas_monopole.txt", status="replace")
            CALL test_forces_time(Nlist, thetas, eps2, tests, file, nt, 10, eval_df, 10)
        CLOSE(file)

        print *, '# quadrupole'
        OPEN(file, file = TRIM(out_dir)//"bh_test_multipoles_thetas_quadrupole.txt", status="replace")
            CALL test_forces_time(Nlist, thetas, eps2, tests, file, nt, 11, eval_df, 10)
        CLOSE(file)

        print *, '# octupole'
        OPEN(file, file = TRIM(out_dir)//"bh_test_multipoles_thetas_octupole.txt", status="replace")
            CALL test_forces_time(Nlist, thetas, eps2, tests, file, nt, 12, eval_df, 10)
        CLOSE(file)
    
    print *, '> DEHNEN METHOD'
        print *, '# monopole'
        OPEN(file, file = TRIM(out_dir)//"dehnen_test_multipoles_thetas_monopole.txt", status="replace")
            CALL test_forces_time(Nlist, thetas, eps2, tests, file, nt, 20, eval_df, 10)
        CLOSE(file)

        print *, '# quadrupole'
        OPEN(file, file = TRIM(out_dir)//"dehnen_test_multipoles_thetas_quadrupole.txt", status="replace")
            CALL test_forces_time(Nlist, thetas, eps2, tests, file, nt, 21, eval_df, 10)
        CLOSE(file)

        print *, '# octupole'
        OPEN(file, file = TRIM(out_dir)//"dehnen_test_multipoles_thetas_octupole.txt", status="replace")
            CALL test_forces_time(Nlist, thetas, eps2, tests, file, nt, 22, eval_df, 10)
        CLOSE(file)

    print *, '> BH GROUPS METHOD'
        print *, '# monopole'
        OPEN(file, file = TRIM(out_dir)//"bhg_test_multipoles_thetas_monopole.txt", status="replace")
            CALL test_forces_time(Nlist, thetas, eps2, tests, file, nt, 30, eval_df, 10)
        CLOSE(file)

        print *, '# quadrupole'
        OPEN(file, file = TRIM(out_dir)//"bhg_test_multipoles_thetas_quadrupole.txt", status="replace")
            CALL test_forces_time(Nlist, thetas, eps2, tests, file, nt, 31, eval_df, 10)
        CLOSE(file)

        print *, '# octupole'
        OPEN(file, file = TRIM(out_dir)//"bhg_test_multipoles_thetas_octupole.txt", status="replace")
            CALL test_forces_time(Nlist, thetas, eps2, tests, file, nt, 32, eval_df, 10)
        CLOSE(file)
END PROGRAM