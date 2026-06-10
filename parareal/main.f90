PROGRAM parareal_test
    USE parareal_mod    
    IMPLICIT NONE

    INTEGER :: N, num_windows, N_itermax, qntt
    REAL(16) :: tf, theta, G, eps, maxerror

    N = 1000
    theta = 1.0d0
    G = 1.0d0
    eps = 0.5d0

    tf = 10.0d0
    num_windows = 30
    N_itermax = 15
    maxerror = 1e-20
    qntt = 16 ! the timestep of the coarse is the fine dt times qntt

    WRITE (*,'(A)') '========================================'
    WRITE (*,'(A12," = ",I10)')      'N', N
    WRITE (*,'(A12," = ",ES12.4)')   'theta', theta
    WRITE (*,'(A12," = ",ES12.4)')   'G', G
    WRITE (*,'(A12," = ",ES12.4)')   'eps', eps
    WRITE (*,*)
    WRITE (*,'(A12," = ",F15.6)')    'tf', tf
    WRITE (*,'(A12," = ",I10)')      'nw', num_windows
    WRITE (*,'(A12," = ",I10)')      'Nitermax', N_itermax
    WRITE (*,*)
    WRITE (*,'(A12," = ",I10)')      'qntt', qntt
    WRITE (*,'(A12," = ",ES12.4)')   'maxerror', maxerror
    WRITE (*,'(A)') '========================================'
    WRITE (*,*)

    CALL test(N, num_windows, N_itermax, tf, theta, G, eps, maxerror, qntt)
CONTAINS

SUBROUTINE generate_initial_values (N, m, qs)
    INTEGER, INTENT(IN) :: N
    REAL(16), INTENT(INOUT) :: m(N), qs(N,3)

    CALL RANDOM_SEED()
    m = 1.0d0 / N

    CALL RANDOM_NUMBER(qs)
    qs = 10.0d0 * (2.0d0 * qs - 1.0d0)
END SUBROUTINE

SUBROUTINE test (N, num_windows, N_itermax, tf, theta, G, eps, maxerror, qntt)
    INTEGER, INTENT(IN) :: N, num_windows, N_itermax, qntt
    REAL(16), INTENT(IN) :: tf, theta, G, eps, maxerror
    REAL(16), ALLOCATABLE :: m(:), q0(:,:), p0(:,:)
    REAL(16), ALLOCATABLE :: q1(:,:,:), p1(:,:,:), u(:,:,:), uk(:,:)
    REAL(16) :: erro, timer_0, timer_1
    INTEGER :: i, j

    ALLOCATE(m(N))
    ALLOCATE(q0(N,3))
    ALLOCATE(p0(N,3))
    p0 = 0.0d0

    CALL generate_initial_values(N, m, q0)
    q0(:,3) = 0.0d0
    
    ALLOCATE(q1(num_windows, N,3))
    ALLOCATE(p1(num_windows, N,3))

    WRITE (*,'(A)') '========================================'
    WRITE (*,'(A)') '> PARAREAL (barnes-hut)'
    WRITE (*,*)

    CALL integrate_parareal(m, q0, p0, tf, num_windows, N_itermax, theta, G, eps, maxerror, qntt, q1, p1)

    WRITE (*,'(A)') '========================================'
    WRITE (*,*)

    WRITE (*,'(A)') '========================================'
    WRITE (*,'(A)') '> PARAREAL (direct)'
    WRITE (*,*)

    CALL integrate_parareal(m, q0, p0, tf, num_windows, N_itermax, -1.0_16, G, eps, maxerror, qntt, q1, p1)

    WRITE (*,'(A)') '========================================'
    WRITE (*,*)

    ! OPEN(23, file = "out/test.txt", status="replace")
    ! DO i = 1, num_windows
    !     WRITE (23, *) i * tf/num_windows, q1(i,:,:)
    ! END DO
    ! CLOSE(23)

    WRITE (*,'(A)') '========================================'
    WRITE (*,'(A)') '> SEQUENTIAL'
    WRITE (*,*)
    ALLOCATE(u(num_windows, 2*N,3))
    ALLOCATE(uk(2*N,3))
    u(1,1:N,:) = q0
    u(1,N+1:,:) = p0
    uk(1:N,:) = q0
    uk(N+1:,:) = p0

    timer_0 = omp_get_wtime()
    DO i = 1, num_windows - 1
        uk = F_til(m, uk, tf/(16*num_windows), 16, theta, G, eps)
        u(i+1,:,:) = uk
    END DO
    timer_1 = omp_get_wtime()

    erro = NORM2(u(num_windows,1:N,:) - q1(num_windows,:,:))**2
    erro = erro + NORM2(u(num_windows,N+1:,:) - p1(num_windows,:,:))**2
    erro = SQRT(erro)/(6*N)

    WRITE(*,'(A15," = ",ES12.4)') '||Par. - Seq.||', erro
    WRITE(*,'(A15," = ",F12.4,A2)') 'Total time', timer_1 - timer_0, 's'

END SUBROUTINE
END PROGRAM