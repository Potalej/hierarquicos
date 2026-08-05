! this program just asserst that the generated using Morton ordering
! is the same as using the usual/intuitive way of constructing octrees
PROGRAM main
    USE morton_octree_mod
    USE octree_mod
    IMPLICIT NONE

    INTEGER, PARAMETER :: pf  = SELECTED_REAL_KIND(15, 307)

    INTEGER :: N
    INTEGER :: L
    REAL(pf), ALLOCATABLE :: masses(:), positions(:,:), forces_morton(:,:)
    REAL(pf), ALLOCATABLE :: forces_intuitive(:,:)
    REAL(pf) :: t0, tf, eps, theta2
    INTEGER :: a, nt
    TYPE(octreetype),       ALLOCATABLE :: intuitive_tree
    TYPE(morton_tree_type), ALLOCATABLE :: morton_tree

    N = 10
    eps = 0.01_pf
    L = 0
    theta2 = 0.1_pf**2
    nt = 10 ! number of threads

    ALLOCATE(masses(N))
    ALLOCATE(positions(3,N))

    masses = 1.0_pf

    positions(:,1)  = [0.45931272_pf, 0.22315240_pf, 0.0_pf]
    positions(:,2)  = [0.40157369_pf, 0.21211072_pf, 0.0_pf]
    positions(:,3)  = [0.70372114_pf, 0.64428738_pf, 0.0_pf]
    positions(:,4)  = [0.20444018_pf, 0.68014233_pf, 0.0_pf]
    positions(:,5)  = [0.73388286_pf, 0.59264995_pf, 0.0_pf]
    positions(:,6)  = [0.68856885_pf, 0.48249009_pf, 0.0_pf]
    positions(:,7)  = [0.48843330_pf, 0.28037004_pf, 0.0_pf]
    positions(:,8)  = [0.15927911_pf, 0.69084037_pf, 0.0_pf]
    positions(:,9)  = [0.30486200_pf, 0.50604824_pf, 0.0_pf]
    positions(:,10) = [0.59676033_pf, 0.80020091_pf, 0.0_pf]

    ! First we generate using the morton module
    ALLOCATE(morton_tree)
    CALL morton_tree % init(N, masses, positions, L, 1.2_pf, .true.)
    ! evaluate the forces
    forces_morton = morton_tree % forces_par(1.0_pf, eps**2, theta2, .true., nt)


    ! Now we generate using the intuitive way
    ALLOCATE(intuitive_tree)
    ALLOCATE(forces_intuitive(N,3))

    CALL intuitive_tree % init(masses, positions(1,:),positions(2,:),positions(3,:))
    CALL intuitive_tree % evaluate_quad()
    
    DO a = 1, N
        forces_intuitive(:,a) = intuitive_tree % forces(a, theta2, 1.0_pf, eps)
    END DO

    print *, "Difference between Morton and intuitive:", NORM2(forces_morton - forces_intuitive)
CONTAINS

FUNCTION compute_forces_direct (m, x, y, z, G, eps) RESULT (forces)
    REAL(pf), INTENT(IN) :: m(:), x(:), y(:), z(:)
    REAL(pf), INTENT(IN) :: G, eps
    INTEGER :: a, b
    REAL(pf) :: dx, dy, dz, dist2, invdist3, f
    REAL(pf) :: forces(3, SIZE(m))

    forces = 0.0_pf

    DO a = 2, SIZE(m)
        DO b = 1, a - 1
            dx = x(b) - x(a)
            dy = y(b) - y(a)
            dz = z(b) - z(a)
            dist2 = dx*dx + dy*dy + dz*dz + eps*eps
            invdist3 = 1.0_pf / SQRT(dist2)
            invdist3 = invdist3 * invdist3 * invdist3
            f = G * m(a) * m(b) * invdist3

            ! forces over 'a'
            forces(1,a) = forces(1,a) + f * dx
            forces(2,a) = forces(2,a) + f * dy
            forces(3,a) = forces(3,a) + f * dz

            ! forces over 'b'
            forces(1,b) = forces(1,b) - f * dx
            forces(2,b) = forces(2,b) - f * dy
            forces(3,b) = forces(3,b) - f * dz
        END DO
    END DO
END FUNCTION

END PROGRAM