MODULE test_forces_mod
    USE octree_mod
    USE omp_lib
    USE initial_values_mod
    USE direct_forces_mod
    IMPLICIT NONE
    PUBLIC
    PRIVATE pf

    INTEGER, PARAMETER :: pf  = SELECTED_REAL_KIND(15, 307)

CONTAINS

SUBROUTINE uniform_list (Nmin, Nmax, Nstep, Nlist)
    INTEGER, INTENT(IN) :: Nmin, Nmax, Nstep
    INTEGER, ALLOCATABLE, INTENT(OUT) :: Nlist(:)
    INTEGER :: qntt, N, i

    qntt = 0
    DO N = Nmin, Nmax, Nstep
        qntt = qntt + 1
    END DO
    ALLOCATE(Nlist(qntt))
    i = 1
    DO N = Nmin, Nmax, Nstep
        Nlist(i) = N
        i = i + 1
    END DO
END SUBROUTINE

SUBROUTINE log_list (Nmin, Nmax, qntt, Nlist)
    INTEGER, INTENT(IN) :: Nmin, Nmax, qntt
    INTEGER, ALLOCATABLE, INTENT(OUT) :: Nlist(:)
    DOUBLE PRECISION :: logmin, logmax, dlog
    INTEGER :: i

    ALLOCATE(Nlist(qntt))

    IF (qntt == 1) THEN
        Nlist(1) = Nmin
        RETURN
    END IF

    logmin = LOG10(DBLE(Nmin))
    logmax = LOG10(DBLE(Nmax))
    dlog   = (logmax - logmin) / DBLE(qntt - 1)

    DO i = 1, qntt
        Nlist(i) = NINT(10.0D0**(logmin + DBLE(i - 1) * dlog))
    END DO

    ! Garante os extremos exatos, sem erro de arredondamento
    Nlist(1)    = Nmin
    Nlist(qntt) = Nmax
END SUBROUTINE

SUBROUTINE cross_product (u, v, w)
    REAL(pf), INTENT(IN) :: u(3), v(3)
    REAL(pf), INTENT(OUT) :: w(3)

    w(1) =  u(2)*v(3)-v(2)*u(3)
    w(2) = -u(1)*v(3)+v(1)*u(3)
    w(3) =  u(1)*v(2)-v(1)*u(2)
END SUBROUTINE

SUBROUTINE evaluate_error_on_accelerations (ms, qs, forces, forces_dir, nt, errors)
    REAL(pf), INTENT(IN)  :: ms(:), qs(:,:)
    REAL(pf), INTENT(OUT) :: errors(:)
    REAL(pf), INTENT(IN)  :: forces(:,:), forces_dir(:,:)
    INTEGER,  INTENT(IN)  :: nt
    REAL(pf) :: accel(3, SIZE(ms)), accel_dir(3, SIZE(ms))
    INTEGER  :: p, N
    REAL(pf) :: a_scale, torque(3)

    N = size(ms)

    ! ERRORS
    ! errors(1,:) -> errors wrt direct forces
    ! errors(2-4,:) -> SUM F_i
    ! errors(5-7,:) -> torque SUM q_i x F_i
    
    a_scale = 0.0_pf
    errors = 0.0_pf

    ! sum of the forces
    errors(2) = SUM(forces(1,:))
    errors(3) = SUM(forces(2,:))
    errors(4) = SUM(forces(3,:))

    IF (nt == 1) THEN
        DO p = 1, N
            ! transforming to accelerations
            accel(:,p) = forces(:,p) / ms(p)
            accel_dir(:,p) = forces_dir(:,p) / ms(p)

            ! scale
            a_scale = a_scale + ms(p) / DOT_PRODUCT(qs(:,p), qs(:,p))

            ! torques
            CALL cross_product(qs(:,p), forces(:,p), torque)
            errors(5) = errors(5) + torque(1)
            errors(6) = errors(6) + torque(2)
            errors(7) = errors(7) + torque(3)
        END DO
    ELSE
        !$OMP PARALLEL DO PRIVATE(p, torque) &
        !$OMP SHARED(ms, forces, forces_dir, accel, accel_dir) &
        !$OMP REDUCTION(+:a_scale,errors) NUM_THREADS(nt)
        DO p = 1, N
            ! transforming to accelerations
            accel(:,p) = forces(:,p) / ms(p)
            accel_dir(:,p) = forces_dir(:,p) / ms(p)

            ! scale
            a_scale = a_scale + ms(p) / DOT_PRODUCT(qs(:,p), qs(:,p))

            ! torques
            CALL cross_product(qs(:,p), forces(:,p), torque)
            errors(5) = errors(5) + torque(1)
            errors(6) = errors(6) + torque(2)
            errors(7) = errors(7) + torque(3)
        END DO
        !$OMP END PARALLEL DO
    ENDIF

    errors(1) = SQRT((NORM2(accel - accel_dir)/a_scale)**2 / N)
END SUBROUTINE

SUBROUTINE test_forces_time (N_list, thetas, eps2, tests, file, nt, method, eval_df, nt_df, gs_par)
    INTEGER,  INTENT(IN) :: N_list(:), tests, file, nt, method
    REAL(pf), INTENT(IN) :: thetas(:), eps2
    LOGICAL,  INTENT(IN) :: eval_df
    INTEGER,  INTENT(IN) :: nt_df
    INTEGER,  INTENT(IN), OPTIONAL :: gs_par

    INTEGER :: N, i_N, i_test, i_theta, p, gs
    REAL(pf), ALLOCATABLE :: m(:), qs(:,:), ps(:,:)       ! state vectors
    REAL(pf), ALLOCATABLE :: forces(:,:), forces_dir(:,:) ! forces

    REAL(pf) :: time_start, time_finish ! timers
    REAL(pf) :: time_pre_init, time_init, time_df, time_forces

    REAL(pf) :: theta, errors(7)
    CLASS(OctreeType), ALLOCATABLE :: tree

    gs = 64
    IF (PRESENT(gs_par)) gs = gs_par

    DO i_N = 1, SIZE(N_list)
        N = N_list(i_N)

        PRINT *, 'N=', N
        ALLOCATE(m(N))
        ALLOCATE(qs(3,N), ps(3,N))
        ALLOCATE(forces(3,N), forces_dir(3,N))

        DO i_test = 1, tests
            CALL generate_initial_values(N, m, qs, ps, &
                7.685706306597839_pf, & ! size_qs
                0.0_pf,  & ! size_ps
                .TRUE.,  & ! normalized
                .FALSE., & ! 2d
                "plummer", & ! mass density profile
                (/1.0_pf, 0.5890486225480862_pf /), & ! profile parameters
                i_test & ! seed
            )

            ! compute forces directly
            IF (eval_df) THEN
                time_start = omp_get_wtime()
                    forces_dir = compute_forces_direct(m, qs(1,:), qs(2,:), qs(3,:), 1.0_pf, eps2, nt_df)
                time_finish = omp_get_wtime()
                time_df = time_finish - time_start
                errors = 0.0_pf
                WRITE(file, *) N, -1.0_pf, 0.0_pf, 0.0_pf, time_df, errors
            ELSE
                forces_dir = 0.0_pf
            ENDIF
            IF (method == 0) CYCLE

            ! starting the octree
            ALLOCATE(tree)
            
            time_start = omp_get_wtime()
                CALL tree % pre_init(m, method, group_size=gs)
            time_finish = omp_get_wtime()
            time_pre_init = time_finish - time_start

            time_start = omp_get_wtime()
                CALL tree % init(qs(1,:), qs(2,:), qs(3,:))
            time_finish = omp_get_wtime()
            time_init = time_finish - time_start

            DO i_theta = 1, SIZE(thetas)
                theta = thetas(i_theta)**2

                SELECT CASE (method)
                    ! classical Barnes-Hut
                    CASE (10, 11, 12)
                        time_start = omp_get_wtime()
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

                    ! Dehnen
                    CASE (20, 21, 22)
                        time_start = omp_get_wtime()
                            CALL tree % dehnen_eval(theta, eps2, 1.0_pf, forces)
                        time_finish = omp_get_wtime()

                    ! Optimized Barnes-Hut
                    CASE (30, 31, 32)
                        time_start = omp_get_wtime()
                            CALL tree % forces_groups(theta, 1.0_pf, eps2, forces, nt)
                        time_finish = omp_get_wtime()

                    CASE DEFAULT
                        STOP "method not identified!"
                END SELECT

                time_forces = time_finish - time_start

                CALL evaluate_error_on_accelerations(m, qs, forces, forces_dir, nt_df, errors)

                WRITE(file, *) N, thetas(i_theta), time_pre_init, time_init, time_forces, errors
            END DO

            DEALLOCATE(tree)
        END DO
        DEALLOCATE(m, qs, ps, forces, forces_dir)
    END DO
END SUBROUTINE

END MODULE