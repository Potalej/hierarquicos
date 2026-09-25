! ************************************************************
!! Initial values
!
!> Objectives
!  This module provides routines to generate random initial
!  values.
!
!> Modified
!  2026.09.24
!
!> Created
!  2026.07.28
!
!> Author
!  oap
!
MODULE initial_values_mod
    IMPLICIT NONE
    PUBLIC generate_initial_values, api
    PRIVATE

    REAL(8), PARAMETER :: PI = 4.0d0 * ATAN(1.0d0)
CONTAINS

SUBROUTINE api (size_qs, size_ps, nm, td, pr, pp, m, qs, ps)
    REAL(8), INTENT(IN) :: size_qs, size_ps
    LOGICAL, INTENT(IN) :: nm ! normalized masses
    LOGICAL, INTENT(IN) :: td ! two dimensional
    CHARACTER(LEN=*), INTENT(IN) :: pr ! mass profile
    REAL(8),          INTENT(IN) :: pp(:) ! profile parameters
    REAL(8), INTENT(INOUT) :: m(:), qs(:,:), ps(:,:)

    CALL generate_initial_values (SIZE(m), m, qs, ps, size_qs, size_ps, nm, td, pr, pp)
END SUBROUTINE

SUBROUTINE generate_initial_values (N, m, qs, ps, size_qs, size_ps, nm, td, pr, pp)
    INTEGER,  INTENT(IN) :: N
    REAL(8), INTENT(IN) :: size_qs, size_ps
    LOGICAL,  INTENT(IN) :: nm ! normalized masses
    LOGICAL,  INTENT(IN) :: td ! two dimensional
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: pr ! mass profile
    REAL(8), INTENT(IN), OPTIONAL         :: pp(:) ! profile parameters
    REAL(8), INTENT(INOUT) :: m(N), qs(3,N), ps(3,N)

    REAL(8) :: tlm(3), com(3)
    INTEGER :: p

    CALL RANDOM_SEED()

    ! if profile not present, generate uniform
    IF (.NOT. PRESENT(pr)) THEN
        CALL generate_initial_values_uniform(N, m, qs, ps, size_qs, size_ps, nm, td)
        RETURN
    ENDIF

    SELECT CASE (TRIM(pr))
        CASE ("uniform")
            CALL generate_initial_values_uniform(N, m, qs, ps, size_qs, size_ps, nm, td)
        
        CASE ("plummer")
            IF (.NOT. PRESENT(pp) .OR. SIZE(pp) .NE. 2) THEN
                STOP "For the Plummer profile, it is necessary to set the parameters G and b."
            ENDIF
            CALL generate_initial_values_plummer(N, m, qs, ps, size_qs, size_ps, nm, td, pp)

        CASE ("homogeneous")
            IF (.NOT. PRESENT(pp) .OR. SIZE(pp) .NE. 3) THEN
                STOP "For the homogeneous sphere profile, it is necessary to set the parameters G, R and rho0."
            ENDIF
            CALL generate_initial_values_homogeneous_sphere(N, m, qs, ps, size_qs, size_ps, &
                                                            nm, td, pp)
        
        CASE ("hernquist")
            IF (.NOT. PRESENT(pp) .OR. SIZE(pp) .NE. 2) THEN
                STOP "For the Hernquist sphere profile, it is necessary to set the parameters G and a"
            ENDIF
            CALL generate_initial_values_hernquist_isotropic(N, m, qs, ps, size_qs, size_ps, &
                                                            nm, td, pp)

        CASE DEFAULT
            STOP 'Unknow density profile: "'//TRIM(pr)//'"'
    END SELECT

    ! total linear momentum -> 0
    tlm(1) = sum(ps(1,:))/N
    tlm(2) = sum(ps(2,:))/N
    tlm(3) = sum(ps(3,:))/N
    DO p = 1, N
        ps(:,p) = ps(:,p) - tlm
    END DO

    ! center of mass to the origin
    com = 0.0d0
    DO p = 1, N
        com = com + m(p) * qs(:,p)
    END DO
    com = com / SUM(m)

    DO p = 1, N
        qs(:,p) = qs(:,p) - com
    END DO
END SUBROUTINE

SUBROUTINE generate_initial_values_uniform (N, m, qs, ps, size_qs, size_ps, nm, two_dim)
! Generates random initial values with uniform distribution (positions
! and momenta). The masses can be 1/N if nm == .TRUE. or uniformly
! distributed with M = 1.
    INTEGER,  INTENT(IN) :: N
    REAL(8), INTENT(IN) :: size_qs, size_ps
    LOGICAL,  INTENT(IN) :: nm, two_dim
    REAL(8), INTENT(INOUT) :: m(N), qs(3,N), ps(3,N)

    CALL random_seed()

    ! normalized masses
    IF (nm) THEN
        m = 1.0d0 / N
    ELSE
        CALL random_number(m)
        DO WHILE (MINVAL(m) == 0.0d0)
            CALL random_number(m)
        END DO
        m = m / SUM(m)
    ENDIF

    CALL random_number(qs)
    qs = size_qs * (2.0d0 * qs - 1.0d0)

    CALL random_number(ps)
    ps = size_ps * (2.0d0 * ps - 1.0d0)

    IF (two_dim) THEN
        qs(3,:) = 0.0d0
        ps(3,:) = 0.0d0
    ENDIF
END SUBROUTINE

SUBROUTINE generate_initial_values_plummer (N, m, qs, ps, size_qs, size_ps, nm, td, pars)
! Generates random initial values with Plummer density and isotropic velocities. 
! The masses can be 1/N if nm == .TRUE. or uniformly distributed with M = 1.
    INTEGER,  INTENT(IN) :: N
    LOGICAL,  INTENT(IN) :: nm, td ! normalized masses, two dimensional
    REAL(8), INTENT(IN) :: pars(2) ! G and plummer parameter
    REAL(8), INTENT(IN) :: size_qs, size_ps
    REAL(8), INTENT(INOUT) :: m(N), qs(3,N), ps(3,N)

    REAL(8) :: X1, X2(N), X3(N)
    REAL(8) :: r(N), phi(N), cos_theta(N), sin_theta(N)
    
    INTEGER  :: i
    REAL(8) :: plummer_potential, v_escape(N)
    REAL(8) :: G, plupar
    REAL(8) :: max_ps(3), big_ps

    REAL(8) :: v, v_phi, v_cos_theta, v_sin_theta
    REAL(8) :: qvn, yvn

    ! parameters
    G = pars(1)
    plupar = pars(2)

    CALL random_seed()

    !> MASSES
    !  m = 1/N
        m = 1.0d0 / N

    !> POSITIONS
    !  Inverting the cumulative distribution function P to get r = P(u), u \sim U[0,1]
        i = 1
        DO WHILE (i <= N)
            ! uniformly distributed values
            CALL random_number(X1)

            ! apply P^-1(u) = r
            r(i) = X1**(1.0d0/3.0d0) * plupar / (SQRT(1.0d0 - X1**(2.0d0/3.0d0)))
            IF (size_qs >= 0 .AND. r(i) > size_qs) CYCLE

            i = i + 1
        END DO

        ! phi angle variables
        CALL random_number(X2)
        phi = 2.0d0 * PI * X2

        ! theta angle variables
        ! if two dimensional, theta = PI
        IF (td) THEN
            cos_theta = 0.0d0
            sin_theta = 1.0d0
        ELSE
            CALL random_number(X3)
            cos_theta = 2.0d0 * X3 - 1.0d0
            sin_theta = SQRT(1.0d0 - cos_theta*cos_theta)
        ENDIF

        ! now evaluate the positions
        qs(1,:) = r * sin_theta * COS(phi)
        qs(2,:) = r * sin_theta * SIN(phi)
        qs(3,:) = r * cos_theta

    !> MOMENTA
    !  By von Neumann rejection
        ! first we evaluate the escape velocities
        DO i = 1, N
            plummer_potential = - G * 1.0d0 / SQRT(r(i)**2 + plupar**2)
            v_escape(i) = SQRT(-2.0d0 * plummer_potential)

            von_neumann: DO WHILE (.TRUE.)
                CALL random_number(qvn)
                CALL random_number(yvn)
                yvn = 0.1d0 * yvn

                IF (yvn < qvn*qvn*(1.0d0 - qvn*qvn)**(3.5d0)) EXIT von_neumann
            END DO von_neumann

            v = qvn * v_escape(i)

            CALL random_number(v_phi)
            v_phi = 2.0d0 * PI * v_phi

            IF (td) THEN
                v_cos_theta = 0.0d0
                v_sin_theta = 1.0d0
            ELSE
                CALL random_number(v_cos_theta)
                v_cos_theta = 2.0d0 * v_cos_theta - 1.0d0
                v_sin_theta = SQRT(1.0d0 - v_cos_theta**2)
            ENDIF

            ps(1,i) = m(i) * v * v_sin_theta * COS(v_phi)
            ps(2,i) = m(i) * v * v_sin_theta * SIN(v_phi)
            ps(3,i) = m(i) * v * v_cos_theta
        END DO
        IF (size_ps >= 0) THEN
            max_ps(1) = MAXVAL(ABS(ps(1,:)))
            max_ps(2) = MAXVAL(ABS(ps(2,:)))
            max_ps(3) = MAXVAL(ABS(ps(3,:)))
            big_ps = MAXVAL(max_ps)

            IF (big_ps > size_ps) ps = ps * size_ps / big_ps
        ENDIF

    !> 2d
    IF (td) THEN
        ps(3,:) = 0.0d0
        qs(3,:) = 0.0d0
    ENDIF
END SUBROUTINE

SUBROUTINE generate_initial_values_homogeneous_sphere (N, m, qs, ps, &
            size_qs, size_ps, nm, td, pars)
! Generates random initial values with homogeneous density and isotropic velocities.
! The masses can be 1/N if nm == .TRUE. or uniformly distributed with M = 1.
    INTEGER,  INTENT(IN) :: N
    LOGICAL,  INTENT(IN) :: nm, td ! normalized masses, two dimensional
    REAL(8), INTENT(IN) :: pars(3) ! G, R and rho0
    REAL(8), INTENT(IN) :: size_qs, size_ps
    REAL(8), INTENT(INOUT) :: m(N), qs(3,N), ps(3,N)

    REAL(8) :: X1, X2(N), X3(N)
    REAL(8) :: r(N), phi(N), cos_theta(N), sin_theta(N)
    
    INTEGER  :: i
    REAL(8) :: potential, v_escape(N)
    REAL(8) :: G, Rpar, rho0par ! parameters of the profile
    REAL(8) :: max_ps(3), big_ps

    ! parameters
    G = pars(1)
    Rpar = pars(2)
    rho0par = pars(3)

    CALL random_seed()

    !> MASSES
    !  m = 1/N
        m = 1.0d0 / N

    !> POSITIONS
    !  Inverting the cumulative distribution function P to get r = P(u), u \sim U[0,1]
        i = 1
        DO WHILE (i <= N)
            ! uniformly distributed values
            CALL random_number(X1)

            ! apply P^-1(u) = r
            r(i) = Rpar * X1 ** (1./3.)
            IF (size_qs >= 0 .AND. r(i) > size_qs) CYCLE

            i = i + 1
        END DO

        ! phi angle variables
        CALL random_number(X2)
        phi = 2.0d0 * PI * X2

        ! theta angle variables
        ! if two dimensional, theta = PI
        IF (td) THEN
            cos_theta = 0.0d0
            sin_theta = 1.0d0
        ELSE
            CALL random_number(X3)
            cos_theta = 2.0d0 * X3 - 1.0d0
            sin_theta = SQRT(1.0d0 - cos_theta*cos_theta)
        ENDIF

        ! now evaluate the positions
        qs(1,:) = r * sin_theta * COS(phi)
        qs(2,:) = r * sin_theta * SIN(phi)
        qs(3,:) = r * cos_theta

    !> MOMENTA
    !  By von Neumann rejection
        ! first we evaluate the escape velocities
        DO i = 1, N
            potential = (2.0d0 * PI * G * rho0par) / 3.0d0
            IF (r(i) < Rpar) THEN
                potential = potential * (r(i)**2 - 3.0d0 * Rpar**2)
            ELSE
                potential = - 2.0d0 * potential * (Rpar**3 / r(i))
            ENDIF
            v_escape(i) = SQRT(-2.0d0 * potential)

            ! TODO
        END DO

    !> 2d
    IF (td) THEN
        ps(3,:) = 0.0d0
        qs(3,:) = 0.0d0
    ENDIF
END SUBROUTINE

SUBROUTINE generate_initial_values_hernquist_isotropic (N, m, qs, ps, &
            size_qs, size_ps, nm, td, pars)
! Generates random initial values with Hernquist density and isotropic velocities.
! The masses can be 1/N if nm == .TRUE. or uniformly distributed with M = 1.
    INTEGER,  INTENT(IN) :: N
    LOGICAL,  INTENT(IN) :: nm, td ! normalized masses, two dimensional
    REAL(8), INTENT(IN) :: pars(3) ! G, R and rho0
    REAL(8), INTENT(IN) :: size_qs, size_ps
    REAL(8), INTENT(INOUT) :: m(N), qs(3,N), ps(3,N)

    REAL(8) :: X1, X2(N), X3(N)
    REAL(8) :: r(N), phi(N), cos_theta(N), sin_theta(N)
    
    INTEGER  :: i
    REAL(8) :: potential, v_escape(N)
    REAL(8) :: G, apar ! parameters of the profile
    REAL(8) :: max_ps(3), big_ps

    ! parameters
    G = pars(1)
    apar = pars(2)

    CALL random_seed()

    !> MASSES
    !  m = 1/N
        m = 1.0d0 / N

    !> POSITIONS
    !  Inverting the cumulative distribution function P to get r = P(u), u \sim U[0,1]
        i = 1
        DO WHILE (i <= N)
            ! uniformly distributed values
            CALL random_number(X1)

            ! apply P^-1(u) = r
            r(i) = apar * SQRT(X1)/(1.0d0 - SQRT(X1))
            IF (size_qs >= 0 .AND. r(i) > size_qs) CYCLE

            i = i + 1
        END DO

        ! phi angle variables
        CALL random_number(X2)
        phi = 2.0d0 * PI * X2

        ! theta angle variables
        ! if two dimensional, theta = PI
        IF (td) THEN
            cos_theta = 0.0d0
            sin_theta = 1.0d0
        ELSE
            CALL random_number(X3)
            cos_theta = 2.0d0 * X3 - 1.0d0
            sin_theta = SQRT(1.0d0 - cos_theta*cos_theta)
        ENDIF

        ! now evaluate the positions
        qs(1,:) = r * sin_theta * COS(phi)
        qs(2,:) = r * sin_theta * SIN(phi)
        qs(3,:) = r * cos_theta

    !> MOMENTA
    !  By von Neumann rejection
        ! first we evaluate the escape velocities
        DO i = 1, N
            potential = - G/ (apar + r(i))
            v_escape(i) = SQRT(-2.0d0 * potential)

            ! TODO
        END DO
        
    !> 2d
    IF (td) THEN
        ps(3,:) = 0.0d0
        qs(3,:) = 0.0d0
    ENDIF
END SUBROUTINE

END MODULE