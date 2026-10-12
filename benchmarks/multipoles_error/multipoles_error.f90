! THE OBJECTIVE HERE IS TO TEST THE ERROR OF THE INTERACTION BETWEEN A PARTICLE AND A NODE
! Ill test the three options of multipole to see to the error varies within the angle parameter 
! and the distance. 
PROGRAM multipoles_error
    USE octree_mod
    USE omp_lib
    USE initial_values_mod
    USE forces_mod
    IMPLICIT NONE

    INTEGER, PARAMETER :: pf  = SELECTED_REAL_KIND(15, 307)
    INTEGER :: N_node
    REAL(pf), ALLOCATABLE :: m(:), qs(:,:), ps(:,:), thetas(:), eps2_list(:)
    REAL(pf) :: min_d, theta, dvec(3), desc, s, q0(3)
    REAL(pf) :: com(3)
    INTEGER :: p, i_theta, num_tests, test, i_eps2
    CLASS(OctreeType), ALLOCATABLE :: tree_mono, tree_quad, tree_octu
    REAL(pf), ALLOCATABLE :: forces_dir(:,:)
    REAL(pf) :: forces_tree(3,3)
    REAL(pf) :: eps2
    INTEGER  :: file_out = 43

    CALL SYSTEM("mkdir -p out")

    N_node = 3000
    num_tests = 20
    ALLOCATE(m(N_node + 1), qs(3, N_node+1), ps(3, N_node+1))
    ALLOCATE(forces_dir(3,N_node+1))

    m(1) = 1.0_pf

    ! given a theta, the distance needed to the sphere be a node is d > s / theta.
    ALLOCATE(thetas(5))
    thetas = [0.2_pf, 0.4_pf, 0.6_pf, 0.8_pf, 1.0_pf]

    ALLOCATE(eps2_list(5))
    eps2_list = [0.0_pf, 0.01_pf**2, 0.05_pf**2, 0.1_pf**2, 0.15_pf**2]
    OPEN(file_out, file="out/multipoles_error.txt", status="replace")

    DO test = 1, num_tests
        ! first, Ill sample a tiny Plummer Sphere that will be the node
        CALL generate_initial_values(N_node, m(2:), qs(:, 2:), ps(:, 2:), &
            6.0_pf, & ! size_qs
            0.0_pf,  & ! size_ps
            .TRUE.,  & ! normalized
            .FALSE., & ! 2d
            "plummer", & ! mass density profile
            (/1.0_pf, 0.58_pf /) & ! profile parameters
        )

        ! now I construct the tree
        ALLOCATE(tree_mono, tree_quad, tree_octu)

        ! for the monopole
        CALL tree_mono % pre_init(m(2:), 10)
        CALL tree_mono % init(qs(1,2:), qs(2,2:), qs(3,2:))

        ! for the quadrupole
        CALL tree_quad % pre_init(m(2:), 11)
        CALL tree_quad % init(qs(1,2:), qs(2,2:), qs(3,2:))

        ! for the octupole
        CALL tree_octu % pre_init(m(2:), 12)
        CALL tree_octu % init(qs(1,2:), qs(2,2:), qs(3,2:))

        qs(:,1) = 1.0_pf
        desc = NORM2(qs(:,1))
        s = SQRT(tree_mono % ns_L2(1))
        
        DO WHILE (s > desc * MINVAL(thetas))
            qs(:,1) = 1.1 * qs(:,1)
            desc = NORM2(qs(:,1))
            print *, desc, s
        END DO
        qs(:,1) = 10.0_pf * qs(:,1)
        q0 = qs(:,1)

        DO i_eps2 = 1, SIZE(eps2_list)
            eps2 = eps2_list(i_eps2)

            DO i_theta = 1, SIZE(thetas)
                theta = thetas(i_theta)
                qs(:,1) = q0

                DO
                    qs(:,1) = 0.9_pf * qs(:,1)
                    desc = NORM2(qs(:,1))
                    print *, "S:", s, desc, desc * theta

                    IF (s > desc * theta) EXIT

                    forces_dir = compute_forces_direct(m, qs(1,:), qs(2,:), qs(3,:), 1.0_pf, eps2, 10)

                    CALL eval_forces_particle_tree(m(1), qs(:,1), tree_mono, 1.0_pf, eps2, 1, forces_tree(:,1))
                    CALL eval_forces_particle_tree(m(1), qs(:,1), tree_quad, 1.0_pf, eps2, 4, forces_tree(:,2))
                    CALL eval_forces_particle_tree(m(1), qs(:,1), tree_octu, 1.0_pf, eps2, 8, forces_tree(:,3))

                    WRITE (file_out, *) 1, theta, SQRT(eps2), desc, NORM2(forces_tree(:,1)-forces_dir(:,1))/NORM2(forces_dir(:,1))
                    WRITE (file_out, *) 4, theta, SQRT(eps2), desc, NORM2(forces_tree(:,2)-forces_dir(:,1))/NORM2(forces_dir(:,1))
                    WRITE (file_out, *) 8, theta, SQRT(eps2), desc, NORM2(forces_tree(:,3)-forces_dir(:,1))/NORM2(forces_dir(:,1))
                END DO

            END DO

        END DO

        DEALLOCATE(tree_mono, tree_quad, tree_octu)

    END DO

    stop 0

    ! and evaluate the forces of 1 over the rest
    ! forces_tree(:,1) = tree % forces(1, theta, 1.0_pf, eps2)

    print *, forces_tree(:,1) - forces_dir(:,1)

CONTAINS

SUBROUTINE eval_forces_particle_tree (m, q, tree, G, eps2, multipole, forces)
    REAL(pf), INTENT(IN) :: m, q(3), G, eps2
    INTEGER,  INTENT(IN) :: multipole
    CLASS(OctreeType), INTENT(IN) :: tree
    REAL(pf), INTENT(OUT) :: forces(3)

    REAL(pf) :: dx, dy, dz, dist2, rinv, f
    REAL(pf) :: invR, invR3, invR5
    REAL(pf) :: Mx, My, Mz
    REAL(pf) :: Qx, Qy, Qz
    REAL(pf) :: mxi2, myi2, mzi2, mxyi, mxzi, myzi

    REAL(pf) :: dx2, dy2, dz2, dxyz, invR2

    REAL(pf) :: invR7, invR9, Ox, Oy, Oz, Oc
    REAL(pf) :: Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz

    forces = 0.0_pf

    dx = tree % ns_qcm_x(1) - q(1)
    dy = tree % ns_qcm_y(1) - q(2)
    dz = tree % ns_qcm_z(1) - q(3)
    dist2 = dx*dx + dy*dy + dz*dz + eps2

    IF (multipole == 1) THEN
        rinv = 1.0_pf / SQRT(dist2)
        rinv = rinv * rinv * rinv
        f = G * m * tree % ns_mass(1) * rinv
        
        forces(1) = f * dx
        forces(2) = f * dy
        forces(3) = f * dz

        RETURN
    ENDIF

    invR = 1.0_pf / SQRT(dist2)
    invR2 = 1.0_pf / dist2
    invR3 = invR * invR2
    invR5 = invR3 * invR * invR
    invR7 = invR5 * invR * invR
    invR9 = invR7 * invR * invR

    Mx = dx * invR3 * tree % ns_mass(1)
    My = dy * invR3 * tree % ns_mass(1)
    Mz = dz * invR3 * tree % ns_mass(1)
    
    mxi2 = tree % ns_quad(1,1)
    myi2 = tree % ns_quad(1,2)
    mzi2 = tree % ns_quad(1,3)
    mxyi = tree % ns_quad(1,4)
    mxzi = tree % ns_quad(1,5)
    myzi = tree % ns_quad(1,6)

    dx2  = 15.0_pf * dx * dx * invR2 - 3.0_pf
    dy2  = 15.0_pf * dy * dy * invR2 - 3.0_pf
    dz2  = 15.0_pf * dz * dz * invR2 - 3.0_pf
    dxyz = 15.0_pf * dx * dy * dz * invR2
    
    Qx = (mxi2 * (dx2 - 6.0_pf) + myi2 * dy2 + mzi2 * dz2) * 0.5_pf * dx
    Qy = (myi2 * (dy2 - 6.0_pf) + mzi2 * dz2 + mxi2 * dx2) * 0.5_pf * dy
    Qz = (mzi2 * (dz2 - 6.0_pf) + mxi2 * dx2 + myi2 * dy2) * 0.5_pf * dz

    Qx = Qx + mxyi * dy * dx2 + mxzi * dz * dx2 + myzi * dxyz
    Qy = Qy + myzi * dz * dy2 + mxyi * dx * dy2 + mxzi * dxyz
    Qz = Qz + mxzi * dx * dz2 + myzi * dy * dz2 + mxyi * dxyz

    ! quadrupole
    IF (multipole == 4) THEN
        forces(1) = G * m * (Mx + invR5 * Qx)
        forces(2) = G * m * (My + invR5 * Qy)
        forces(3) = G * m * (Mz + invR5 * Qz)
    
    ! octupole
    ELSE
        Ox = - 1.5_pf * invR5 * tree % ns_oct(1, 12)
        Oy = - 1.5_pf * invR5 * tree % ns_oct(1, 13)
        Oz = - 1.5_pf * invR5 * tree % ns_oct(1, 14)

        Oc = 7.5_pf * invR7 * (dx * tree % ns_oct(1, 12) + &
                               dy * tree % ns_oct(1, 13) + &
                               dz * tree % ns_oct(1, 14))
        Ox = Ox + dx * Oc
        Oy = Oy + dy * Oc
        Oz = Oz + dz * Oc

        dx2 = dx*dx
        dy2 = dy*dy
        dz2 = dz*dz
        Ox = Ox + 7.5_pf * invR7 * (dx2 * tree % ns_oct(1, 1) + &
                                    dy2 * tree % ns_oct(1, 2) + &
                                    dz2 * tree % ns_oct(1, 3) + &
                                    2.0_pf * dx*dy * tree % ns_oct(1, 4) + &
                                    2.0_pf * dx*dz * tree % ns_oct(1, 7) + &
                                    2.0_pf * dy*dz * tree % ns_oct(1, 10))
        Oy = Oy + 7.5_pf * invR7 * (dx2 * tree % ns_oct(1, 4) + &
                                    dy2 * tree % ns_oct(1, 5) + &
                                    dz2 * tree % ns_oct(1, 6) + &
                                    2.0_pf * dx*dy * tree % ns_oct(1, 2) + &
                                    2.0_pf * dx*dz * tree % ns_oct(1, 10) + &
                                    2.0_pf * dy*dz * tree % ns_oct(1, 8))
        Oz = Oz + 7.5_pf * invR7 * (dx2 * tree % ns_oct(1, 7) + &
                                    dy2 * tree % ns_oct(1, 8) + &
                                    dz2 * tree % ns_oct(1, 9) + &
                                    2.0_pf * dx*dy * tree % ns_oct(1, 10) + &
                                    2.0_pf * dx*dz * tree % ns_oct(1, 3) + &
                                    2.0_pf * dy*dz * tree % ns_oct(1, 6))

        Oc = -17.5_pf * invR9 * ( &
            dx2 * dx * tree % ns_oct(1,1) + &
            dy2 * dy * tree % ns_oct(1,5) + &
            dz2 * dz * tree % ns_oct(1,9) + &
            3.0_pf * dx2 * dy * tree % ns_oct(1,4) + &
            3.0_pf * dx2 * dz * tree % ns_oct(1,7) + &
            3.0_pf * dy2 * dx * tree % ns_oct(1,2) + &
            3.0_pf * dy2 * dz * tree % ns_oct(1,8) + &
            3.0_pf * dz2 * dx * tree % ns_oct(1,3) + &
            3.0_pf * dz2 * dy * tree % ns_oct(1,6) + &
            6.0_pf * dx * dy * dz * tree % ns_oct(1,10))

        Ox = Ox + dx * Oc
        Oy = Oy + dy * Oc
        Oz = Oz + dz * Oc

        forces(1) = forces(1) + G * m * (Mx + invR5 * Qx + Ox)
        forces(2) = forces(2) + G * m * (My + invR5 * Qy + Oy)
        forces(3) = forces(3) + G * m * (Mz + invR5 * Qz + Oz)
    ENDIF
END SUBROUTINE

END PROGRAM