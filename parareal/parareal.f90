MODULE parareal_mod
USE api_mod
USE OMP_LIB
IMPLICIT NONE
PUBLIC
CONTAINS

FUNCTION G_til (m, u0, dt, qntt, theta, G, eps)
    REAL(16), INTENT(IN) :: m(:), u0(:,:), dt, theta, G, eps
    INTEGER, INTENT(IN) :: qntt
    REAL(16) :: G_til(SIZE(u0,1),3)
    
    REAL(16) :: forces(INT(SIZE(u0,1)/2),3)
    REAL(16) :: q(INT(SIZE(u0,1)/2),3), p(INT(SIZE(u0,1)/2),3)
    INTEGER :: N, i
    
    N = SIZE(q,1)
    q = u0(1:N,:)
    p = u0(N+1:,:)

    IF (theta == -1) THEN
        DO i=1, qntt
            forces = compute_forces_direct(m, q(:,1), q(:,2), q(:,3), G, eps)
            p = p + dt * forces
            q = q + dt * p
        END DO
    ELSE
        DO i=1, qntt
            forces = compute_forces(m, q(:,1), q(:,2), q(:,3), theta, G, eps, .TRUE.)
            p = p + dt * forces
            q = q + dt * p
        END DO
    ENDIF

    G_til(1:N,:) = q
    G_til(N+1:,:) = p
END FUNCTION

FUNCTION F_til (m, u0, dt, qntt, theta, G, eps)
    REAL(16), INTENT(IN) :: m(:), u0(:,:), dt, theta, G, eps
    REAL(16) :: F_til(SIZE(u0,1),3)
    REAL(16) :: forces(INT(SIZE(u0,1)/2),3)
    REAL(16) :: q(INT(SIZE(u0,1)/2),3), p(INT(SIZE(u0,1)/2),3)
    REAL(16) :: p_half(INT(SIZE(u0,1)/2),3)
    INTEGER :: qntt, i, N
    
    N = SIZE(q,1)
    q = u0(1:N,:)
    p = u0(N+1:,:)

    forces = compute_forces_direct(m, q(:,1), q(:,2), q(:,3), G, eps)
    DO i=1, qntt
        ! verlet
        p_half = p + 0.5d0 * dt * forces
        q = q + dt * p_half
        forces = compute_forces_direct(m, q(:,1), q(:,2), q(:,3), G, eps)
        p = p_half + 0.5d0 * dt * forces
    END DO

    F_til(1:N,:) = q
    F_til(N+1:,:) = p
END FUNCTION

FUNCTION F_til_PAR (m, u0, dt, qntt, theta, G, eps)
    REAL(16), INTENT(IN) :: m(:), u0(:,:), dt, theta, G, eps
    REAL(16) :: F_til_PAR(SIZE(u0,1),3)
    REAL(16) :: forces(INT(SIZE(u0,1)/2),3)
    REAL(16) :: q(INT(SIZE(u0,1)/2),3), p(INT(SIZE(u0,1)/2),3)
    REAL(16) :: p_half(INT(SIZE(u0,1)/2),3)
    INTEGER :: qntt, i, N
    
    N = SIZE(q,1)
    q = u0(1:N,:)
    p = u0(N+1:,:)

    forces = compute_forces_direct_par(m, q(:,1), q(:,2), q(:,3), G, eps)
    DO i=1, qntt
        ! verlet
        p_half = p + 0.5d0 * dt * forces
        q = q + dt * p_half
        forces = compute_forces_direct_par(m, q(:,1), q(:,2), q(:,3), G, eps)
        p = p_half + 0.5d0 * dt * forces
    END DO

    F_til_PAR(1:N,:) = q
    F_til_PAR(N+1:,:) = p
END FUNCTION

SUBROUTINE integrate_parareal (m, q0, p0, tf, num_windows, N_itermax, theta, G, eps, maxerror, qntt, q1, p1)
    REAL(16), INTENT(IN) :: m(:), q0(:,:), p0(:,:), theta, G, eps
    REAL(16), INTENT(IN) :: tf, maxerror
    INTEGER, INTENT(IN) :: num_windows, N_itermax, qntt
    REAL(16), INTENT(OUT) :: q1(:,:,:), p1(:,:,:)

    REAL(16) :: window_size, dt_fino, dt_grosseiro, error
    INTEGER :: i, k, N, iterations
    REAL(16), ALLOCATABLE :: U(:,:,:,:), U_til(:,:,:), u0(:,:)
    REAL(16) :: timer_0, timer_1, timer_total

    window_size = tf / num_windows
    dt_grosseiro = window_size
    dt_fino = dt_grosseiro / qntt

    N = SIZE(q0,1)
    ALLOCATE(u0(2*N,3))
    u0(1:N,  :) = q0
    u0(N+1:, :) = p0
    
    ALLOCATE(U(N_itermax, num_windows, 2*N, 3))
    U = 0.0d0

    ALLOCATE(U_til(num_windows, 2*N, 3))

    U(1,1,:,:) = u0
    DO i = 1, num_windows - 1
        U(1,i+1,:,:) = G_til(m, U(1,i,:,:), dt_grosseiro, 1, theta, G, eps)
    END DO

    iterations = 1
    timer_total = 0
    DO k = 1, N_itermax - 1
        iterations = iterations + 1
        timer_0 = omp_get_wtime()

        U(k+1,1,:,:) = u0
        U_til(1,:,:) = u0

        !$OMP PARALLEL PRIVATE(i)
        !$OMP DO
        DO i = 1, num_windows - 1
            U_til(i+1,:,:) = F_til(m, U(k,i,:,:), dt_fino, qntt, theta, G, eps)
        END DO
        !$OMP END DO
        !$OMP END PARALLEL        

        DO i = 1, num_windows - 1
            U(k+1,i+1,:,:) = G_til(m, U(k+1,i,:,:), dt_grosseiro, 1, theta, G, eps) + U_til(i+1,:,:)
            U(k+1,i+1,:,:) = U(k+1,i+1,:,:) - G_til(m, U(k,i,:,:), dt_grosseiro, 1, theta, G, eps)
        END DO

        timer_1 = omp_get_wtime()
        timer_total = timer_total + timer_1 - timer_0

        error = NORM2(U(k+1, num_windows,:,:) - U(k, num_windows,:,:))/(6*N)
        ! error = NORM2(U(k+1, num_windows,:,:) - U_til(num_windows,:,:))/(6*N)

        WRITE(*,'(A6,I4,A9,ES12.4,A8,ES12.4)') "Iter.:", iterations, " / Error:", error, " / Time:", timer_1 - timer_0
        IF (error <= maxerror) EXIT
    END DO

    WRITE(*,*)
    WRITE(*,'(A17,ES12.4)') 'Max Error Accep.', maxerror
    WRITE(*,'(A17,ES12.4)') 'Final error', error
    WRITE(*,'(A17,F12.4,A2)') 'Total time', timer_total, 's'

    q1 = U(iterations,:,1:N,:)
    p1 = U(iterations,:,N+1:,:)

    DEALLOCATE(U, U_til)
END SUBROUTINE

END MODULE