PROGRAM integration
    USE octree_mod
    USE omp_lib
    USE initial_values_mod
    USE forces_mod
    IMPLICIT NONE

    INTEGER, PARAMETER :: pf  = SELECTED_REAL_KIND(15, 307)
    INTEGER :: seed_size
    INTEGER, ALLOCATABLE :: seed(:)
    INTEGER :: N, nt, multipole
    REAL(pf) :: t0, t1, tf, dt, eps2, G, E, E0, theta2, p0(3)
    REAL(pf), ALLOCATABLE :: ms(:), qs(:,:), ps(:,:), qs0(:,:), ps0(:,:)

    CALL RANDOM_SEED(size=seed_size)
    ALLOCATE(seed(seed_size))
    seed = 12345
    CALL RANDOM_SEED(put=seed)

    N = 100
    G = 1.0_pf
    eps2 = 0.05_pf ** 2
    dt = 0.0078125_pf
    tf = 1.0_pf
    nt = 1
    theta2 = 0.5_pf**2
    multipole = 1

    ALLOCATE(ms(N), qs(3,N), ps(3,N), qs0(3,N), ps0(3,N))

    CALL generate_initial_values(N, ms, qs0, ps0, &
        6.5238301325205335_pf, & ! size_qs
        -1.0_pf,  & ! size_ps
        .TRUE.,  & ! normalized
        .FALSE., & ! 2d
        "plummer", & ! mass density profile
        (/1.0_pf, 0.59_pf /) & ! profile parameters
    )

    CALL total_energy(ms, qs0, ps0, G, 0.0_pf, E0)
    p0(1) = SUM(ps0(1,:))
    p0(2) = SUM(ps0(2,:))
    p0(3) = SUM(ps0(3,:))

    qs = qs0
    ps = ps0
    t0 = omp_get_wtime()
    CALL test_integration_verlet_direct(ms, qs, ps, G, eps2, dt, tf, nt)
    t1 = omp_get_wtime()
    CALL total_energy(ms, qs, ps, G, 0.0_pf, E)
    PRINT *, 'E:', E, E0, ABS(E - E0)
    PRINT *, 'P:', NORM2((/sum(ps(1,:))-p0(1), sum(ps(2,:))-p0(2), sum(ps(3,:))-p0(3)/))
    PRINT *, 'tempo:', t1 - t0

    print * , ''
    print * , '-----------'
    print * , ''

    qs = qs0
    ps = ps0
    t0 = omp_get_wtime()
    CALL test_integration_verlet_bh(ms, qs, ps, G, eps2, dt, tf, theta2, multipole, nt)
    t1 = omp_get_wtime()
    CALL total_energy(ms, qs, ps, G, eps2, E)
    PRINT *, 'E:', E, E0, ABS(E - E0)
    PRINT *, 'P:', NORM2((/sum(ps(1,:))-p0(1), sum(ps(2,:))-p0(2), sum(ps(3,:))-p0(3)/))
    PRINT *, 'tempo:', t1 - t0

    print * , ''
    print * , '-----------'
    print * , ''

    qs = qs0
    ps = ps0
    t0 = omp_get_wtime()
    CALL test_integration_verlet_dehnen(ms, qs, ps, G, eps2, dt, tf, theta2, multipole)
    t1 = omp_get_wtime()
    CALL total_energy(ms, qs, ps, G, eps2, E)
    PRINT *, 'E:', E, E0, ABS(E - E0)
    PRINT *, 'P:', NORM2((/sum(ps(1,:))-p0(1), sum(ps(2,:))-p0(2), sum(ps(3,:))-p0(3)/))
    PRINT *, 'tempo:', t1 - t0

CONTAINS

SUBROUTINE test_integration_verlet_direct (ms, qs, ps, G, eps2, dt, tf, nt)
    REAL(pf), INTENT(INOUT) :: qs(:,:), ps(:,:)
    REAL(pf), INTENT(IN) :: ms(:), G, dt, tf, eps2
    INTEGER :: nt

    REAL(pf) :: m
    REAL(pf) :: fs(3,SIZE(ms))
    REAL(pf) :: t

    m = ms(1)

    t = 0.0_pf
    fs = compute_forces_direct(ms, qs(1,:), qs(2,:), qs(3,:), G, eps2, nt)
    DO WHILE (t < tf)
        ps = ps + 0.5_pf * dt * fs
        qs = qs + dt * ps / m

        fs = compute_forces_direct(ms, qs(1,:), qs(2,:), qs(3,:), G, eps2, nt)

        ps = ps + 0.5_pf * dt * fs

        t = t + ABS(dt)
    END DO
END SUBROUTINE

SUBROUTINE test_integration_verlet_bh (ms, qs, ps, G, eps2, dt, tf, theta2, mult, nt)
    REAL(pf), INTENT(INOUT) :: qs(:,:), ps(:,:)
    REAL(pf), INTENT(IN) :: ms(:), G, dt, tf, eps2, theta2
    INTEGER,  INTENT(IN) :: nt, mult

    INTEGER  :: p, N
    REAL(pf) :: m
    REAL(pf) :: fs(3,SIZE(ms))
    REAL(pf) :: t
    INTEGER  :: method = 10

    CLASS(OctreeType), ALLOCATABLE :: tree

    N = SIZE(ms)
    m = ms(1)

    t = 0.0_pf

    ALLOCATE(tree)
    IF (mult == 4) method = 11
    IF (mult == 8) method = 12
    CALL tree % pre_init(ms, method)

    CALL tree % init(qs(1,:), qs(2,:), qs(3,:))

    IF (nt == 1) THEN
        DO p = 1, N
            fs(:,p) = tree % forces(p, theta2, G, eps2)
        END DO
    ELSE
        !$OMP PARALLEL DO SHARED(fs) PRIVATE(p) NUM_THREADS(nt) &
        !$OMP SCHEDULE(DYNAMIC)
        DO p = 1, N
            fs(:,p) = tree % forces(p, theta2, G, eps2)
        END DO
        !$OMP END PARALLEL DO
    END IF


    DO WHILE (t < tf)
        ps = ps + 0.5_pf * dt * fs
        qs = qs + dt * ps / m

        CALL tree % init(qs(1,:), qs(2,:), qs(3,:))

        IF (nt == 1) THEN
            DO p = 1, N
                fs(:,p) = tree % forces(p, theta2, G, eps2)
            END DO
        ELSE
            !$OMP PARALLEL DO SHARED(fs) PRIVATE(p) NUM_THREADS(nt) &
            !$OMP SCHEDULE(DYNAMIC)
            DO p = 1, N
                fs(:,p) = tree % forces(p, theta2, G, eps2)
            END DO
            !$OMP END PARALLEL DO
        END IF

        ps = ps + 0.5_pf * dt * fs

        t = t + ABS(dt)
    END DO
END SUBROUTINE

SUBROUTINE test_integration_verlet_dehnen (ms, qs, ps, G, eps2, dt, tf, theta2, mult)
    REAL(pf), INTENT(INOUT) :: qs(:,:), ps(:,:)
    REAL(pf), INTENT(IN) :: ms(:), G, dt, tf, eps2, theta2
    INTEGER,  INTENT(IN) :: mult

    INTEGER  :: N
    REAL(pf) :: m
    REAL(pf) :: fs(3,SIZE(ms))
    REAL(pf) :: t
    INTEGER  :: method = 20

    CLASS(OctreeType), ALLOCATABLE :: tree

    N = SIZE(ms)
    m = ms(1)
    t = 0.0_pf

    ALLOCATE(tree)
    IF (mult == 4) method = 21
    CALL tree % pre_init(ms, method)

    CALL tree % init(qs(1,:), qs(2,:), qs(3,:))
    CALL tree % dehnen_eval(theta2, eps2, G, fs)

    DO WHILE (t < tf)
        ps = ps + 0.5_pf * dt * fs
        qs = qs + dt * ps / m

        CALL tree % init(qs(1,:), qs(2,:), qs(3,:))
        CALL tree % dehnen_eval(theta2, eps2, G, fs)

        ps = ps + 0.5_pf * dt * fs

        t = t + ABS(dt)
    END DO
END SUBROUTINE

SUBROUTINE total_energy (ms, qs, ps, G, eps2, E)
    REAL(pf), INTENT(IN) :: ms(:), qs(:,:), ps(:,:), G, eps2
    REAL(pf) :: E, dist, T, V
    INTEGER  :: a, b

    T = DOT_PRODUCT(ps(:,1), ps(:,1))/(2.0_pf * ms(1))
    V = 0.0_pf
    
    DO a = 2, SIZE(ms)
        T = T + DOT_PRODUCT(ps(:,a), ps(:,a))/(2.0_pf * ms(a))
        
        DO b = 1, a - 1
            dist = SQRT((qs(1,b) - qs(1,a))**2 + (qs(2,b) - qs(2,a))**2 + (qs(3,b) - qs(3,a))**2 + eps2)
            V = V - G * ms(a) * ms(b) / dist
        END DO
    END DO

    ! print *, T, V, E
    E = T + V
END SUBROUTINE

END PROGRAM