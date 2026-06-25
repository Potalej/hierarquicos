! this program just asserst that the generated using Morton ordering
! is the same as using the usual/intuitive way of constructing octrees
PROGRAM main
    USE morton_octree_mod
    USE octree_mod
    IMPLICIT NONE

    INTEGER :: N
    INTEGER :: L
    REAL(8), ALLOCATABLE :: masses(:), positions(:,:), forces_morton(:,:)
    REAL(8), ALLOCATABLE :: forces_intuitive(:,:)
    REAL(8) :: t0, tf, eps, theta2
    INTEGER :: a
    TYPE(octreetype),       ALLOCATABLE :: intuitive_tree
    TYPE(morton_tree_type), ALLOCATABLE :: morton_tree

    N = 10
    eps = 0.01d0
    L = 0
    theta2 = 0.1d0**2

    ALLOCATE(masses(N))
    ALLOCATE(positions(N,3))

    masses = 1.0d0

    positions(1,:)  = [0.45931272d0, 0.22315240d0, 0.0d0]
    positions(2,:)  = [0.40157369d0, 0.21211072d0, 0.0d0]
    positions(3,:)  = [0.70372114d0, 0.64428738d0, 0.0d0]
    positions(4,:)  = [0.20444018d0, 0.68014233d0, 0.0d0]
    positions(5,:)  = [0.73388286d0, 0.59264995d0, 0.0d0]
    positions(6,:)  = [0.68856885d0, 0.48249009d0, 0.0d0]
    positions(7,:)  = [0.48843330d0, 0.28037004d0, 0.0d0]
    positions(8,:)  = [0.15927911d0, 0.69084037d0, 0.0d0]
    positions(9,:)  = [0.30486200d0, 0.50604824d0, 0.0d0]
    positions(10,:) = [0.59676033d0, 0.80020091d0, 0.0d0]

    ! First we generate using the morton module
    ALLOCATE(morton_tree)
    CALL morton_tree % init(N, masses, positions, L, 1.2d0, .true.)
    ! evaluate the forces
    forces_morton = morton_tree % forces(1.0d0, eps**2, theta2, .true.)


    ! Now we generate using the intuitive way
    ALLOCATE(intuitive_tree)
    ALLOCATE(forces_intuitive(N,3))

    CALL intuitive_tree % init(masses, positions(:,1),positions(:,2),positions(:,3))
    CALL intuitive_tree % evaluate_quad()
    
    DO a = 1, N
        forces_intuitive(a,:) = intuitive_tree % forces(a, theta2, 1.0d0, eps)
    END DO

    print *, "Difference between Morton and intuitive:", NORM2(forces_morton - forces_intuitive)
CONTAINS

FUNCTION compute_forces_direct (m, x, y, z, G, eps) RESULT (forces)
    REAL(8), INTENT(IN) :: m(:), x(:), y(:), z(:)
    REAL(8), INTENT(IN) :: G, eps
    INTEGER :: a, b
    REAL(8) :: dx, dy, dz, dist2, invdist3, f
    REAL(8) :: forces(SIZE(m), 3)

    forces = 0.0d0

    DO a = 2, SIZE(m)
        DO b = 1, a - 1
            dx = x(b) - x(a)
            dy = y(b) - y(a)
            dz = z(b) - z(a)
            dist2 = dx*dx + dy*dy + dz*dz + eps*eps
            invdist3 = 1.0d0 / SQRT(dist2)
            invdist3 = invdist3 * invdist3 * invdist3
            f = G * m(a) * m(b) * invdist3

            ! forces over 'a'
            forces(a,1) = forces(a,1) + f * dx
            forces(a,2) = forces(a,2) + f * dy
            forces(a,3) = forces(a,3) + f * dz

            ! forces over 'b'
            forces(b,1) = forces(b,1) - f * dx
            forces(b,2) = forces(b,2) - f * dy
            forces(b,3) = forces(b,3) - f * dz
        END DO
    END DO
END FUNCTION

END PROGRAM