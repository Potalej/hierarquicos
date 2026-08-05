! ************************************************************
!! Octree for Barnes-Hut using Morton (Z) ordering
!
!> Objectives
!  This module generates an octree based in the positions by
!  normalizing all the system into the cube [0,2^L-1]^3 and 
!  generating a Morton key for each body. The forces are
!  evaluated using the Barnes-Hut method with a quadrupole
!  expansion.
!
!> Modified
!  2026.08.03
!
!> Created
!  2026.08.03
!
!> Author
!  oap
!
MODULE morton_octree_mod
    USE omp_lib
    IMPLICIT NONE
    PRIVATE
    PUBLIC morton_tree_type

    INTEGER, PARAMETER :: pf  = SELECTED_REAL_KIND(15, 307)

    TYPE :: node_type
        INTEGER(8) :: prefix ! hash prefix
        INTEGER    :: level  ! 0, ..., L, where 0 is the root
        REAL(pf)    :: size2  ! squared length of the cube
        
        INTEGER :: first ! first particle index (ordered)
        INTEGER :: last  ! last particle index (ordered)

        INTEGER :: parent           ! node parent
        INTEGER :: nchild, child(8) ! children
        INTEGER :: next             ! next node index

        REAL(pf) :: mass   ! mass
        REAL(pf) :: qcm(3) ! center of mass
        
        REAL(pf) :: quad(6) ! quadrupole expansion
    END TYPE

    TYPE :: morton_tree_type
        INTEGER :: N ! number of bodies
        TYPE(node_type), ALLOCATABLE :: nodes(:)
        INTEGER :: nnodes
        
        INTEGER, ALLOCATABLE :: perm(:)
        REAL(pf), ALLOCATABLE :: masses(:), positions(:,:)
        REAL(pf), ALLOCATABLE :: perm_masses(:), perm_positions(:,:)
        INTEGER, ALLOCATABLE :: level_begin(:), level_end(:)
        REAL(pf) :: system_size_amplificator

        INTEGER :: depth
        LOGICAL :: use_quadrupole

        REAL(pf) :: time_generate_tree, time_eval_quadrupole

        CONTAINS
            PROCEDURE :: init => init_morton_tree_type
            PROCEDURE :: forces_over_p => forces_individual_morton_tree_type
            PROCEDURE :: forces_seq
            PROCEDURE :: forces_par
    END TYPE

CONTAINS

! MORTON_TREE_TYPE PROCEDURES
! init the class, generate the tree and evaluate the octree if its the case
SUBROUTINE init_morton_tree_type (self, N, ms, qs, L, ssa, use_quadrupole)
    CLASS(morton_tree_type), INTENT(INOUT) :: self
    INTEGER, INTENT(IN) :: N
    REAL(pf), INTENT(IN) :: ms(:), qs(:,:)
    INTEGER, INTENT(IN) :: L
    REAL(pf), INTENT(IN) :: ssa ! system size amplificator
    LOGICAL, INTENT(IN) :: use_quadrupole
    REAL(pf) :: t0

    self % N = N
    ALLOCATE(self % masses(N))
    ALLOCATE(self % positions(3,N))
    self % masses = ms
    self % positions = qs
    self % system_size_amplificator = ssa
    self % depth = L

    ALLOCATE(self % perm(N))
    ALLOCATE(self % perm_masses(N))
    ALLOCATE(self % perm_positions(3,N))

    self % use_quadrupole = use_quadrupole

    ! generate the tree
    t0 = omp_get_wtime()
    CALL generate_nodes(N, ms, qs, self % depth, self % nodes, self % nnodes, &
                        self % perm, self % level_begin, self % level_end, &
                        self % perm_masses, self % perm_positions, &
                        self % system_size_amplificator)
    self % time_generate_tree = omp_get_wtime() - t0
    
    ! evaluate the quadrupoles if wanted
    IF (self % use_quadrupole) THEN
        t0 = omp_get_wtime()
        CALL evaluate_quadrupoles(self % nodes, self % nnodes, self % perm_masses, self % perm_positions, &
                                self % depth, self % level_begin, self % level_end)
        self % time_eval_quadrupole = omp_get_wtime() - t0
    ENDIF
END SUBROUTINE

! it evaluates the individual forces over a particle
FUNCTION forces_individual_morton_tree_type (self, p, G, eps2, theta2, mult_par) RESULT(forces)
    CLASS(morton_tree_type), INTENT(INOUT) :: self
    INTEGER, INTENT(IN) :: p
    REAL(pf), INTENT(IN) :: G, eps2, theta2
    LOGICAL, INTENT(IN), OPTIONAL :: mult_par

    REAL(pf) :: forces(3)
    LOGICAL :: mult

    mult = self % use_quadrupole
    IF (PRESENT(mult_par)) mult = mult_par

    forces = barnes_hut(p, self % nodes, self % perm_masses, self % perm_positions, &
                        theta2, G, eps2, mult)
END FUNCTION

! it evaluates the forces over all particles (sequential)
FUNCTION forces_seq (self, G, eps2, theta2, mult_par) RESULT(forces)
    CLASS(morton_tree_type), INTENT(INOUT) :: self
    REAL(pf), INTENT(IN) :: G, eps2, theta2
    LOGICAL, INTENT(IN), OPTIONAL :: mult_par

    REAL(pf) :: forces(3,self%N)
    LOGICAL :: mult
    INTEGER :: p

    mult = self % use_quadrupole
    IF (PRESENT(mult_par)) mult = mult_par

    DO p = 1, self % N
        forces(:,self % perm(p)) = barnes_hut(p, self % nodes, self % perm_masses, self % perm_positions, &
                            theta2, G, eps2, mult)
    END DO
END FUNCTION

! it evaluates the forces over all particles (parallel)
FUNCTION forces_par (self, G, eps2, theta2, mult_par, num_threads) RESULT(forces)
    CLASS(morton_tree_type), INTENT(INOUT) :: self
    REAL(pf), INTENT(IN) :: G, eps2, theta2
    LOGICAL, INTENT(IN), OPTIONAL :: mult_par
    INTEGER, INTENT(IN), OPTIONAL :: num_threads

    REAL(pf) :: forces(3,self%N)
    LOGICAL :: mult
    INTEGER :: p

    mult = self % use_quadrupole
    IF (PRESENT(mult_par)) mult = mult_par

    IF (PRESENT(num_threads)) THEN
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p) &
        !$OMP NUM_THREADS(num_threads)
        DO p = 1, self % N
            forces(:, self % perm(p)) = barnes_hut(p, self % nodes, self % perm_masses, self % perm_positions, &
                                theta2, G, eps2, mult)
        END DO
        !$OMP END PARALLEL DO
    ELSE
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p)
        DO p = 1, self % N
            forces(:, self % perm(p)) = barnes_hut(p, self % nodes, self % perm_masses, self % perm_positions, &
                                theta2, G, eps2, mult)
        END DO
        !$OMP END PARALLEL DO
    ENDIF    

END FUNCTION

! GERAL SUBROUTINES
! for a given (x,y,z) and depth L, it gives a Morton key
FUNCTION morton_key (q3, L) RESULT(key)
    REAL(pf), INTENT(IN) :: q3(3)
    INTEGER, INTENT(IN) :: L

    REAL(pf)    :: x, y, z
    INTEGER(8) :: ix, iy, iz
    INTEGER(8) :: xb, yb, zb
    INTEGER(8) :: scale, b, oct
    INTEGER(8) :: key

    x = q3(1)
    y = q3(2)
    z = q3(3)

    ! scale = 2**L
    scale = SHIFTL(1_8, L)
    ix = MIN(FLOOR(x * scale), scale-1)
    iy = MIN(FLOOR(y * scale), scale-1)
    iz = MIN(FLOOR(z * scale), scale-1)

    key = 0_8
    DO b = L-1, 0, -1
        xb = IAND(SHIFTR(ix, b), 1_8) ! (ix >> b) & 1
        yb = IAND(SHIFTR(iy, b), 1_8) ! (iy >> b) & 1
        zb = IAND(SHIFTR(iz, b), 1_8) ! (iz >> b) & 1

        oct = ISHFT(xb,2) + ISHFT(yb,1) + zb
        key = ISHFT(key,3) + oct
    END DO
END FUNCTION

! for a given positions vector qs, it gives the minimum depth to
! achieve a tree with only one particle for leaf (brute force)
SUBROUTINE choose_depth (N, qs, L, keys)
    INTEGER, INTENT(IN) :: N
    REAL(pf), INTENT(IN) :: qs(3,N)
    INTEGER,    INTENT(INOUT) :: L       ! depth
    INTEGER(8), INTENT(INOUT) :: keys(N) ! Morton keys with depth L

    INTEGER :: p
    LOGICAL :: not_unique, found

    L = 1
    found = .FALSE.

    search: DO WHILE (.NOT. found)
        DO p = 1, N
            keys(p) = morton_key(qs(:,p), L)

            IF (p > 1) THEN
                not_unique = ANY(keys(1:p-1) == keys(p))
                IF (not_unique) THEN
                    L = L + 1
                    CYCLE search
                ENDIF
            END IF
        END DO

        found = .TRUE.
    END DO search
END SUBROUTINE

! for a given key, a level and a depth L, it gives the prefix
! of the key
FUNCTION prefix (key, level, L)
    INTEGER(8), INTENT(IN) :: key
    INTEGER,    INTENT(IN) :: level, L
    INTEGER(8) :: prefix

    prefix = SHIFTR(key, 3*(L-level))
END FUNCTION

! for given information, it creates a node (node_type)
FUNCTION create_node (level, prefix, first, last, mass, qcm, size) RESULT(node)
    INTEGER,    INTENT(IN) :: level
    INTEGER(8), INTENT(IN) :: prefix
    INTEGER,    INTENT(IN) :: first, last
    REAL(pf),   INTENT(IN) :: mass, qcm(3), size
    TYPE(node_type) :: node

    node % level = level
    node % prefix = prefix
    node % first = first
    node % last = last
    node % mass = mass
    node % qcm = qcm
    node % size2 = size**2
    node % parent = -1
    node % nchild = 0
    node % child = -1
    node % quad = 0.0_pf
END FUNCTION

! for a given positions vector q(N,3) it gives the size of the minimum
! cube that contains all the particles and its center
FUNCTION system_size_info (N, q) RESULT(infos)
    INTEGER, INTENT(IN) :: N
    REAL(pf), INTENT(IN) :: q(3,N)
    REAL(pf) :: x(N), y(N), z(N)
    REAL(pf) :: xmin, xmax, ymin, ymax, zmin, zmax
    REAL(pf) :: infos(4)

    x = q(1,:)
    y = q(2,:)
    z = q(3,:)

    xmin = MINVAL(x)
    xmax = MAXVAL(x)
    ymin = MINVAL(y)
    ymax = MAXVAL(y)
    zmin = MINVAL(z)
    zmax = MAXVAL(z)

    infos(1) = 0.5_pf * (xmin + xmax)
    infos(2) = 0.5_pf * (ymin + ymax)
    infos(3) = 0.5_pf * (zmin + zmax)
    infos(4) = MAXVAL((/ xmax - xmin, ymax - ymin, zmax - zmin /))
END FUNCTION

! for a given node, get the next node
RECURSIVE SUBROUTINE build_next (idx, sibling, nodes)
    INTEGER, INTENT(IN) :: idx
    INTEGER, INTENT(IN) :: sibling

    TYPE(node_type), INTENT(INOUT) :: nodes(:)

    INTEGER :: i, c

    nodes(idx) % next = sibling

    ! if its empty
    IF (nodes(idx) % nchild == 0) RETURN

    ! for each child (except the last), add the next as other child
    DO i = 1, nodes(idx) % nchild - 1
        c = nodes(idx) % child(i)
        CALL build_next(c, nodes(idx) % child(i+1), nodes)
    END DO

    ! for the last child, the next is the next node
    c = nodes(idx) % child(nodes(idx) % nchild)
    CALL build_next(c, sibling, nodes)
END SUBROUTINE

SUBROUTINE generate_nodes (N, masses, positions, & ! state vectors (IN)
                        L_par,  & ! depth parameter (INOUT)
                        nodes,  & ! nodes vector    (OUT)
                        nnodes, & ! number of nodes (OUT)
                        perm,   & ! permutation vector (INOUT)
                        level_begin, & ! for which index each level begins
                        level_end,   & ! for which index each level ends
                        ms, qs, & ! permuted masses and positions
                        system_size_amplificator & ! normalization constant
                        )
    INTEGER, INTENT(IN) :: N
    REAL(pf), INTENT(IN) :: masses(N), positions(3,N)
    TYPE(node_type), ALLOCATABLE, INTENT(OUT) :: nodes(:)
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: level_begin(:), level_end(:)
    INTEGER, INTENT(INOUT) :: L_par
    INTEGER, INTENT(OUT)   :: nnodes
    INTEGER, INTENT(INOUT) :: perm(N)
    REAL(pf), INTENT(INOUT) :: ms(N), qs(3,N)
    REAL(pf), INTENT(IN)    :: system_size_amplificator

    REAL(pf) :: ssi(4) ! system size info
    REAL(pf) :: ssa ! system size amplified
    REAL(pf) :: qcm(3), mass ! center of mass and mass
    INTEGER :: i, p, level, parent ! auxiliar
    REAL(pf) :: positions_norm(3,N) ! normalized positions in [0,1]^3
    INTEGER :: L ! depth
    INTEGER(8), ALLOCATABLE :: prefixes(:,:) ! prefixes of the keys
    INTEGER(8) :: keys(N), keys_tmp(N) ! Morton keys
    INTEGER(8) :: current, pref ! current prefix
    INTEGER :: child_first, first, last ! first and last particles in a level

    ! get infos about the size of the system to normalize it to [0,1]^3
    ssi = system_size_info(N, positions)
    ssa = system_size_amplificator * ssi(4)

    ! center of mass and total mass of the system
    qcm  = 0.0_pf
    mass = 0.0_pf
    DO p = 1, N
        qcm = qcm + masses(p) * positions(:,p)
        mass = mass + masses(p)
    END DO
    qcm = qcm / mass

    ! normalized coordinates in [0,1]^3
    DO p = 1, N
        positions_norm(:,p) = positions(:,p) - qcm
    END DO
    positions_norm = positions_norm / ssa
    
    ! update the system size info
    ! ssi(1:3) = (ssi(1:3) - qcm)/(mass * ssa)
    ssi = system_size_info(N, positions_norm)
    ! translate the center of the system to (1/2,1/2,1/2)
    positions_norm(1,:) = positions_norm(1,:) + 0.5_pf - ssi(1)
    positions_norm(2,:) = positions_norm(2,:) + 0.5_pf - ssi(2)
    positions_norm(3,:) = positions_norm(3,:) + 0.5_pf - ssi(3)

    ! choose the depth or use the informed deepth
    L = L_par
    IF (L == 0) THEN
        CALL choose_depth(N, positions_norm, L, keys)
        L_par = L
    ELSE
        DO p = 1, N
            keys(p) = morton_key(positions_norm(:,p), L)
        END DO
    ENDIF

    ! init and end index of each level
    IF (ALLOCATED(level_begin)) DEALLOCATE(level_begin)
    IF (ALLOCATED(level_end))   DEALLOCATE(level_end)
    ALLOCATE(level_begin(L+1))
    ALLOCATE(level_end(L+1))

    ! sort the keys
    CALL argsort(keys, perm)

    ! permutation of the masses and positions vector and get the prefixes
    ALLOCATE(prefixes(N,L+1))
    DO p = 1, N
        keys_tmp(p) = keys(perm(p))
        ms(p) = masses(perm(p))
        qs(:,p) = positions(:,perm(p))

        DO level = 0, L
            prefixes(p,level+1) = prefix(keys_tmp(p), level, L)
        END DO
    END DO
    keys = keys_tmp

    ! create the nodes
    IF (ALLOCATED(nodes)) DEALLOCATE(nodes)
    ALLOCATE(nodes((L+1)*N))

    ! for the level zero (root), we add every particle
    nnodes = 1
    nodes(1) = create_node(0, prefixes(1,1), 1, N, mass, qcm, ssa)
    level_begin(1) = 1
    level_end(1)   = 1

    ! create the other levels
    level_loop: DO level = 1, L
        level_begin(level + 1) = nnodes + 1

        parent_loop: DO parent = level_begin(level), level_end(level)
            first = nodes(parent) % first
            last  = nodes(parent) % last
            current = -1_8

            particles_loop: DO i = first, last
                pref = prefixes(i, level + 1)

                ! first particle
                ! update the state variables
                IF (current == -1_8) THEN
                    current = pref
                    child_first = i
                    mass = ms(i)
                    qcm = mass * qs(:,i)

                ! if its a different node
                ELSE IF (pref .NE. current) THEN
                    nnodes = nnodes + 1
                    qcm = qcm / mass

                    ! create the child node
                    nodes(nnodes) = create_node(&
                        level, current, child_first, i-1, mass, qcm, &
                        ssa*(2.0_pf**(-level)) &
                        )
                    nodes(nnodes) % parent = parent
                    
                    ! update the parent
                    nodes(parent) % nchild = nodes(parent) % nchild + 1
                    nodes(parent) % child(nodes(parent) % nchild) = nnodes

                    ! update the state variables
                    current = pref
                    child_first = i
                    mass = ms(i)
                    qcm = mass * qs(:,i)

                ! if its the same node
                ELSE
                    mass = mass + ms(i)
                    qcm = qcm + ms(i) * qs(:,i)
                ENDIF
            END DO particles_loop

            ! create the child node
            nnodes = nnodes + 1
            qcm = qcm / mass
            nodes(nnodes) = create_node(&
                level, current, child_first, last, mass, qcm, &
                ssa*(2.0_pf**(-level)) &
            )
            nodes(nnodes) % parent = parent

            ! update the parent
            nodes(parent) % nchild = nodes(parent) % nchild + 1
            nodes(parent) % child(nodes(parent) % nchild) = nnodes
        END DO parent_loop

        level_end(level + 1) = nnodes
    END DO level_loop

    DEALLOCATE(prefixes)

    ! build the nexts
    CALL build_next(1, 0, nodes)
END SUBROUTINE

SUBROUTINE argsort (keys, perm)
    INTEGER(8), INTENT(IN)  :: keys(:)
    INTEGER,    INTENT(OUT) :: perm(:)

    INTEGER :: i
    INTEGER, ALLOCATABLE :: tmp(:)

    DO i = 1, SIZE(keys)
        perm(i) = i
    END DO

    ALLOCATE(tmp(SIZE(keys)))
    
    CALL mergesort(keys, perm, tmp, 1, SIZE(keys))
    
    DEALLOCATE(tmp)
END SUBROUTINE

RECURSIVE SUBROUTINE mergesort (keys, perm, tmp, left, right)
    INTEGER(8), INTENT(IN)    :: keys(:)
    INTEGER,    INTENT(INOUT) :: perm(:)
    INTEGER,    INTENT(INOUT) :: tmp(:)
    INTEGER,    INTENT(IN)    :: left, right
    INTEGER :: mid, i, j, k

    IF (left >= right) RETURN

    mid = (left + right)/2

    CALL mergesort(keys, perm, tmp, left, mid)    ! sort left piece of array
    CALL mergesort(keys, perm, tmp, mid+1, right) ! sort right piece of array

    i = left
    j = mid+1
    k = left
    DO WHILE (i <= mid .AND. j <= right)
        ! in case of a tie, get the left
        IF (keys(perm(i)) <= keys(perm(j))) THEN
            tmp(k) = perm(i)
            i = i + 1
        ELSE
            tmp(k) = perm(j)
            j = j + 1
        ENDIF
        k = k + 1
    END DO

    DO WHILE (i <= mid)
        tmp(k) = perm(i)
        i = i + 1
        k = k + 1
    END DO
    DO WHILE (j <= right)
        tmp(k) = perm(j)
        j = j + 1
        k = k + 1
    END DO
    
    perm(left:right) = tmp(left:right)
END SUBROUTINE

SUBROUTINE evaluate_quadrupoles (nodes, nnodes, ms, qs, max_level, level_begin, level_end)
    TYPE(node_type), INTENT(INOUT) :: nodes(:)
    INTEGER, INTENT(IN) :: nnodes  ! number of nodes
    INTEGER, INTENT(IN) :: max_level, level_begin(:), level_end(:) ! level information
    REAL(pf), INTENT(IN) :: ms(:), qs(:,:) ! permuted state vectors
    
    INTEGER :: node_idx, child_idx, p, level
    REAL(pf) :: pm, px, py, pz, dxi, dyi, dzi

    INTEGER :: child_pos

    DO level = max_level, 0, -1
        DO node_idx = level_begin(level+1), level_end(level+1)
            nodes(node_idx) % quad = 0.0_pf
            IF (nodes(node_idx) % nchild <= 1) THEN
                DO p = nodes(node_idx) % first, nodes(node_idx) % last
                    pm = ms(p)
                    px = qs(1,p)
                    py = qs(2,p)
                    pz = qs(3,p)

                    dxi = px - nodes(node_idx) % qcm(1)
                    dyi = py - nodes(node_idx) % qcm(2)
                    dzi = pz - nodes(node_idx) % qcm(3)

                    nodes(node_idx) % quad(1) = nodes(node_idx) % quad(1) + pm * dxi**2
                    nodes(node_idx) % quad(2) = nodes(node_idx) % quad(2) + pm * dyi**2
                    nodes(node_idx) % quad(3) = nodes(node_idx) % quad(3) + pm * dzi**2
                    nodes(node_idx) % quad(4) = nodes(node_idx) % quad(4) + pm * dxi * dyi
                    nodes(node_idx) % quad(5) = nodes(node_idx) % quad(5) + pm * dxi * dzi
                    nodes(node_idx) % quad(6) = nodes(node_idx) % quad(6) + pm * dyi * dzi
                END DO
            ELSE
                DO child_pos = 1, nodes(node_idx)%nchild

                    child_idx = nodes(node_idx)%child(child_pos)

                    pm = nodes(child_idx)%mass
                    dxi = nodes(child_idx) % qcm(1) - nodes(node_idx) % qcm(1)
                    dyi = nodes(child_idx) % qcm(2) - nodes(node_idx) % qcm(2)
                    dzi = nodes(child_idx) % qcm(3) - nodes(node_idx) % qcm(3)

                    nodes(node_idx) % quad(1) = nodes(node_idx) % quad(1) + &
                        nodes(child_idx) % quad(1) + pm * dxi * dxi

                    nodes(node_idx) % quad(2) = nodes(node_idx) % quad(2) + &
                        nodes(child_idx) % quad(2) + pm * dyi * dyi

                    nodes(node_idx) % quad(3) = nodes(node_idx) % quad(3) + &
                        nodes(child_idx) % quad(3) + pm * dzi * dzi

                    nodes(node_idx) % quad(4) = nodes(node_idx) % quad(4) + &
                        nodes(child_idx) % quad(4) + pm * dxi * dyi

                    nodes(node_idx) % quad(5) = nodes(node_idx) % quad(5) + &
                        nodes(child_idx) % quad(5) + pm * dxi * dzi

                    nodes(node_idx) % quad(6) = nodes(node_idx) % quad(6) + &
                        nodes(child_idx) % quad(6) + pm * dyi * dzi
                END DO
            ENDIF
        END DO
    END DO
END SUBROUTINE

FUNCTION barnes_hut (p, & ! particle index (sorted)
    nodes, & ! list of nodes
    ms, qs, & ! state vectors (sorted)
    theta2, & ! angle parameter squared
    G, eps2, & ! gravitational constant and softening squared
    mult_par & ! use of multipoles
    ) RESULT(forces)
    TYPE(node_type), INTENT(IN) :: nodes(:)
    INTEGER, INTENT(IN) :: p
    REAL(pf), INTENT(IN) :: ms(:), qs(:,:)
    REAL(pf), INTENT(IN) :: theta2, G, eps2         
    LOGICAL, INTENT(IN), OPTIONAL :: mult_par

    REAL(pf) :: forces(3)
    LOGICAL :: mult ! multipole use
    LOGICAL :: only_one_particle
    REAL(pf) :: mp, qp(3) ! particle state
    INTEGER :: idx, b
    REAL(pf) :: dx, dy, dz, dist2, invdist3, f ! forces
    ! multipole variables
    REAL(pf) :: invR, invR2, invR3, invR5 
    REAL(pf) :: mx, my, mz
    REAL(pf) :: mxi2, myi2, mzi2, mxyi, mxzi, myzi
    REAL(pf) :: dx2, dy2, dz2, dxyz
    REAL(pf) :: Qx, Qy, Qz

    mult = .FALSE.
    IF (PRESENT(mult_par)) mult = mult_par

    mp = ms(p)
    qp = qs(:,p)

    forces = 0.0_pf
    idx = 1

    DO WHILE (idx .NE. 0)
        only_one_particle = (nodes(idx) % first == nodes(idx) % last)

        ! if is there only one particle and its p, cycle
        IF (only_one_particle .AND. nodes(idx) % last == 0) THEN
            idx = nodes(idx) % next
            CYCLE
        ENDIF

        ! if is there no childs, go to the particles inside
        IF (nodes(idx) % nchild == 0) THEN
            DO b = nodes(idx) % first, nodes(idx) % last
                IF (p == b) CYCLE

                dx = qs(1,b) - qp(1)
                dy = qs(2,b) - qp(2)
                dz = qs(3,b) - qp(3)
                dist2 = dx*dx + dy*dy + dz*dz + eps2
                invdist3 = 1.0_pf / SQRT(dist2)
                invdist3 = invdist3 * invdist3 * invdist3
                f = G * mp * ms(b) * invdist3
                
                forces(1) = forces(1) + f * dx
                forces(2) = forces(2) + f * dy
                forces(3) = forces(3) + f * dz
            END DO

            idx = nodes(idx) % next
            CYCLE
        ENDIF

        ! if the particle index is in the node, go to the first child
        IF (nodes(idx) % first <= p .AND. p <= nodes(idx) % last) THEN
            idx = nodes(idx) % child(1)
            CYCLE
        ENDIF

        ! in other cases, or is there only (different) particle, or the BH
        ! criteria is valid or we need to go to the first child
        dx = nodes(idx) % qcm(1) - qp(1)
        dy = nodes(idx) % qcm(2) - qp(2)
        dz = nodes(idx) % qcm(3) - qp(3)
        dist2 = dx*dx + dy*dy + dz*dz

        IF (only_one_particle .OR. nodes(idx) % size2 <= theta2 * dist2) THEN
            IF (nodes(idx) % mass > 0) THEN
                dist2 = dist2 + eps2 ! softening

                ! if not multipole or is there only one particle
                IF (.NOT. mult .OR. only_one_particle) THEN
                    invdist3 = 1.0_pf / SQRT(dist2)
                    invdist3 = invdist3 * invdist3 * invdist3
                    f = G * mp * nodes(idx) % mass * invdist3

                    forces(1) = forces(1) + f * dx
                    forces(2) = forces(2) + f * dy
                    forces(3) = forces(3) + f * dz
                
                ! apply multipole
                ELSE
                    invR  = 1.0_pf / SQRT(dist2)
                    invR2 = 1.0_pf / dist2
                    invR3 = invR * invR2
                    invR5 = invR3 * invR * invR

                    mx = dx * invR3 * nodes(idx) % mass
                    my = dy * invR3 * nodes(idx) % mass
                    mz = dz * invR3 * nodes(idx) % mass

                    mxi2 = nodes(idx) % quad(1)
                    myi2 = nodes(idx) % quad(2)
                    mzi2 = nodes(idx) % quad(3)
                    mxyi = nodes(idx) % quad(4)
                    mxzi = nodes(idx) % quad(5)
                    myzi = nodes(idx) % quad(6)

                    dx2 = 15.0_pf * dx * dx * invR2 - 3.0_pf
                    dy2 = 15.0_pf * dy * dy * invR2 - 3.0_pf
                    dz2 = 15.0_pf * dz * dz * invR2 - 3.0_pf
                    dxyz = 15.0_pf * dx * dy * dz * invR2

                    Qx = (mxi2 * (dx2 - 6.0_pf) + myi2 * dy2 + mzi2 * dz2) * 0.5_pf * dx
                    Qy = (myi2 * (dy2 - 6.0_pf) + mzi2 * dz2 + mxi2 * dx2) * 0.5_pf * dy
                    Qz = (mzi2 * (dz2 - 6.0_pf) + mxi2 * dx2 + myi2 * dy2) * 0.5_pf * dz

                    Qx = Qx + mxyi * dy * dx2 + mxzi * dz * dx2 + myzi * dxyz
                    Qy = Qy + myzi * dz * dy2 + mxyi * dx * dy2 + mxzi * dxyz
                    Qz = Qz + mxzi * dx * dz2 + myzi * dy * dz2 + mxyi * dxyz

                    forces(1) = forces(1) + G * mp * (mx + invR5 * Qx)
                    forces(2) = forces(2) + G * mp * (my + invR5 * Qy)
                    forces(3) = forces(3) + G * mp * (mz + invR5 * Qz)
                ENDIF
            ENDIF

            idx = nodes(idx) % next
            CYCLE
        
        ! go to the first child
        ELSE
            idx = nodes(idx) % child(1)
            CYCLE
        ENDIF
    END DO
END FUNCTION

END MODULE