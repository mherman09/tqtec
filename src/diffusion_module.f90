!--------------------------------------------------------------------------------------------------!
! Module Diffusion                                                                                 !
!                                                                                                  !
! Authors:                                                                                         !
!     - Matt Herman (current Modern Fortran version, i.e., what you are looking at right now!)     !
!     - Rachel Piotraschke (original version; PSU MS thesis)                                       !
!     - Kevin Furlong (original version; supervisor)                                               !
!                                                                                                  !
! This module contains variables and subroutines numerically solving the production-diffusion      !
! differential equation. The derivation of the spherical finite difference approximation and the   !
! original Matlab program HeAge are from Piotraschke (2012) (MS thesis).                           !
!                                                                                                  !
! The 3-D production-diffusion equation in a Cartesian coordinate system can be written as         !
!                                                                                                  !
!       du             d2u      d2u      d2u                                                       !
!      ----  =  D * [ ----  +  ----  +  ---- ]  +  P                                               !
!       dt             dx2      dy2      dz2                                                       !
!                                                                                                  !
! where                                                                                            !
!     u: function that will undergo diffusion (temperature, concentration, topography, etc.)       !
!     t: time                                                                                      !
!     D: diffusivity                                                                               !
!     x, y, z: position                                                                            !
!     P: rate of production                                                                        !
!                                                                                                  !
! For spherical or cylindrical coordinate systems, the right hand side uses the corresponding      !
! Laplacian operator.                                                                              !
!                                                                                                  !
! Assuming diffusion has symmetry so that it only occurs in 1-D, the production-diffusion equation !
! in a Cartesian coordinate system becomes (Carslaw and Jaeger, 1959)                              !
!                                                                                                  !
!       du             d2u                                                                         !
!      ----  =  D * [ ---- ]  +  P                                                                 !
!       dt             dx2                                                                         !
!                                                                                                  !
! Or in a spherical coordinate system:                                                             !
!                                                                                                  !
!       du            1       d           du                                                       !
!      ----  =  D * ----- * ---- [ r^2 * ---- ]  +  P                                              !
!       dt           r^2     dr           dr                                                       !
!                                                                                                  !
! where                                                                                            !
!     r: radius                                                                                    !
!                                                                                                  !
! For 2-D diffusion, the production-diffusion equation in a Cartesian coordinate system becomes    !
!                                                                                                  !
!       du             d2u      d2u                                                                !
!      ----  =  D * [ ----  +  ----  ]  +  P                                                       !
!       dt             dx2      dy2                                                                !
!                                                                                                  !
! Or in a cylindrical coordinate system with different r and z diffusivities:                      !
!                                                                                                  !
!       du              1      d         du                d2u                                     !
!      ----  =  Dr * [ --- * ---- ( r * ---- ) ] + Dz * [ ---- ]  +  P                             !
!       dt              r     dr         dr                dz2                                     !
!                                                                                                  !
! References:                                                                                      !
! Carslaw, H.S., Jaeger, J.C. (1959). Conduction of Heat in Solids. Second Edition. Oxford         !
!     University Press. 310 pp.                                                                    !
! Piotraschke, R. (2012). Thermal and Geologic Constraints on the Cretaceous-to-Neogene Tectonic   !
!     Development of the Klamath Mountains, Northern California. M.Sc. Thesis, Pennsylvania State  !
!     University.                                                                                  !
!--------------------------------------------------------------------------------------------------!


module diffusion

    ! Spherical finite difference parameters
    integer :: nnodes_sphere
    double precision, allocatable :: radius_shell(:)
    double precision, allocatable :: volume_shell(:)
    double precision :: volume_sphere

    ! Subroutines
    PUBLIC :: init_spherical_node_geometry
    PUBLIC :: diffusion_step_dirichlet
    PRIVATE :: tridiag


contains



    subroutine init_spherical_node_geometry(r,dr,n)
    !----
    ! Create spherical node geometry for finite difference procedure with two extra nodes at ends
    ! for boundary condition specifications
    !
    !   Center of sphere
    !          |
    !          |           dr                              \ r
    !          V        -------                             |
    !       o  X  *     *     *     *     *     *     *     *     o
    ! BC node                                               |     BC node
    !                                                      /
    !
    ! Inputs (arguments):
    !   r:   radius of sphere
    !   dr:  spatial step size
    !   n:   number of spatial nodes, including 2 BC nodes
    !
    ! Outputs:
    !   radius_shell:     distance from center of sphere to each node [m]
    !   volume_shell:     volume of shell corresponding to each node  [m^3]
    !   volume_sphere:    volume of entire sphere                     [m^3]
    !----

    implicit none


    ! Arguments
    double precision :: r
    double precision :: dr
    integer :: n

    ! Local variables
    integer :: i


    ! Allocate memory to shell radius/volume arrays
    nnodes_sphere = n
    if (.not.allocated(radius_shell)) then
        allocate(radius_shell(nnodes_sphere))
    endif
    if (.not.allocated(volume_shell)) then
        allocate(volume_shell(nnodes_sphere))
    endif

    ! Distance from center of sphere at each spatial node
    do i = 1,nnodes_sphere
        radius_shell(i) = (dble(i)-1.5d0)*dr
    enddo

    ! Volume of each spherical shell
    volume_shell(1) = 0.0d0
    volume_shell(2) = 4.0d0/3.0d0*3.14159265d0 * radius_shell(2)**3
    do i = 3,nnodes_sphere
        volume_shell(i) = 4.0d0/3.0d0*3.14159265d0*(radius_shell(i)**3-radius_shell(i-1)**3)
    enddo

    ! Sphere volume
    volume_sphere = 4.0d0/3.0d0*3.14159265d0 * r**3

    return
    end subroutine



!--------------------------------------------------------------------------------------------------!



    subroutine diffusion_step_dirichlet(u,n,s,beta)
    !----
    ! Solve the finite difference approximation to the diffusion equation for the value of the
    ! function u(i,j) at the next time step u(i,j+1), given fixed values at boundary nodes
    ! (Dirichlet BCs):
    !     1/k * [u(i,j+1) - u(i,j)] =
    !         alpha/h^2 * [ beta*(u(i-1,j+1) - 2*u(i,j+1) + u(i+1,j+1))
    !                           + (1-beta)*(u(i-1,j) - 2*u(i,j) + u(i+1,j) ]
    !
    !     u: function undergoing diffusion (temperature, concentration, topography, etc.)
    !     n: number of spatial nodes in finite difference grid
    !     k: time step size
    !     h: space step size
    !     i: current spatial node
    !     j: current time node
    !     alpha: diffusivity
    !     beta: implicitness weight
    !
    ! This turns into the following system of equations to solve for time step (j+1):
    !
    !     (1+2*s*beta)*u(i,j+1) - s*beta*(u(i-1,j+1)+u(i+1,j+1))
    !         = (1-2*s*(1-beta))*u(i,j) + s*(1-beta)*(u(i-1,j)+u(i+1,j))
    !
    !     s: diffusion number (diffusivity * k/h^2)
    !----

    implicit none

    ! Arguments
    integer :: n
    double precision :: u(n)
    double precision :: s
    double precision :: beta

    ! Local variables
    integer :: i
    double precision :: a(n,3)
    double precision :: b(n,1)


    ! Load the tridiagonal matrix for the LHS of the system of equations
    do i = 1,n
        a(i,1) = -s*beta
        a(i,2) = 2.0d0*s*beta + 1.0d0
        a(i,3) = -s*beta
    enddo

    ! Boundary condition nodes are fixed values (Dirichlet BCs)
    a(1,1) = 0.0d0
    a(1,2) = 1.0d0
    a(1,3) = 0.0d0
    a(n,1) = 0.0d0
    a(n,2) = 1.0d0
    a(n,3) = 0.0d0

    ! Load the RHS vector for the system of equations
    do i = 2,n-1
        b(i,1) = s*(1.0d0-beta)*u(i-1) + (1.0d0-(2.0d0*s*(1.0d0-beta)))*u(i) + s*(1.0d0-beta)*u(i+1)
    enddo

    ! Boundary condition nodes are fixed values (Dirichlet BCs)
    b(1,1) = u(1)
    b(n,1) = u(n)

    ! Solve for the function u at the next time step
    call tridiag(a(:,1),a(:,2),a(:,3),b(:,1),u,n)


    return
    end subroutine




!--------------------------------------------------------------------------------------------------!


#ifdef COMPILE_WITH_SUPERLU

    subroutine diffusion_2d_cylinder(u,nnodes_r,nnodes_z,dr,dz,diffusivity_r,diffusivity_z,dt,beta)
    !----
    ! Solve the finite difference approximation to the cylindrical diffusion equation for the value
    ! of the function u(i,j,k) at the next time step u(i,j,k+1), given fixed values at boundary
    ! nodes (Dirichlet BCs):
    !
    !     1/k * [u(i,j,k+1) - u(i,j,k)] =
    !             Dr/hr^2 * [ beta*(u(i-1,j,k+1) - 2*u(i,j,k+1) + u(i+1,j,k+1))
    !                           + (1-beta)*(u(i-1,j,k) - 2*u(i,j,k) + u(i+1,j,k) ]
    !             Dr/(2*r*hr) * [ beta*(u(i+1,j,k+1) - u(i-1,j,k+1))
    !                                + (1-beta)*(u(i+1,j,k) - u(i-1,j,k) ]
    !             Dz/hz^2 * [ beta*(u(i,j-1,k+1) - 2*u(i,j,k+1) + u(i,j+1,k+1))
    !                           + (1-beta)*(u(i,j-1,k) - 2*u(i,j,k) + u(i,j+1,k) ]
    !
    !     u: function undergoing diffusion (temperature, concentration, topography, etc.)
    !     i: current radial node (r-dimension)
    !     j: current axial node (z-dimension)
    !     k: time step size
    !     r: radial position
    !     hr: space step size in radial direction
    !     hz: space step size in axial direction
    !     Dr: diffusivity in the radial direction
    !     Dz: diffusivity in the axial direction
    !     beta: implicitness weight
    !
    ! This turns into the following system of equations to solve for time step (k+1):
    !
    !     (1 + 2*sz*beta + 2*sr1*beta) * u(i,   j  , k+1)
    !          + (sr2*beta - sr1*beta) * u(i-1, j  , k+1)
    !          - (sr2*beta + sr1*beta) * u(i+1, j  , k+1)
    !                        - sz*beta * u(i,   j-1, k+1)
    !                        - sz*beta * u(i,   j+1, k+1) = 
    !                                 (1 - 2*sz*beta - 2*sr1*beta) * u(i,j,k) 
    !                                     + sr1*(1-beta) * (u(i-1.j,k) + u(i+1,j,k))
    !                                     +  sz*(1-beta) * (u(i,j-1,k) + u(i,j+1,k))
    !                                     + sr2*(1-beta) * (u(i+1,j,k) - u(i-1,j,k))
    !
    !     sz: axial diffusion number (Dz * k/hz^2)
    !     sr1: radial diffusion number (Dr * k/hr^2)
    !     sr2: a different radial diffusion constant (Dr * k/(2*r*hr))
    !
    ! The system of equations (Ax = b) does not form a trivial tridiagonal A matrix like in 1-D,
    ! but A has many zeros so we can use a sparse matrix solver.
    !----


    implicit none

    ! Arguments
    integer :: nnodes_r
    integer :: nnodes_z
    double precision :: u(nnodes_r,nnodes_z)
    double precision :: diffusivity_r
    double precision :: diffusivity_z
    double precision :: dr
    double precision :: dz
    double precision :: dt
    double precision :: beta

    ! Local variables
    double precision :: r
    double precision :: sz
    double precision :: sr1
    double precision :: sr2
    integer :: i                  ! index for A
    integer :: n                  ! number of rows/columns in A
    integer :: nnz                ! number of non-zero entries in A
    integer :: irow               ! A and b row index 
    integer :: icol               ! A column index 
    integer :: ir                 ! radial node index
    integer :: iz                 ! axial node index
    double precision, allocatable :: a(:)
    double precision, allocatable :: b(:)
    integer, allocatable :: row_indices(:)
    integer, allocatable :: column_pointers(:)


    sz  = diffusivity_z * dt / dz**2
    sr1 = diffusivity_r * dt / dr**2

    n = nnodes_r*nnodes_z
    nnz = 5*nnodes_r*nnodes_z ! NEEDS TO BE LARGER FOR BC NODES

    ! Define sparse matrix arrays
    allocate(a(nnz))
    allocate(b(n))
    allocate(row_indices(nnz))
    allocate(column_pointers(n+1))
    a = 0.0d0
    b = 0.0d0
    row_indices = 0
    column_pointers = 0

    ! Set A matrix sparse structure and load values
    i = 1
    do icol = 1,n ! Sparse matrix storage is column-based

        column_pointers(icol) = i

        do iz = 1,nnodes_z ! Loop through all rows in the column
            do ir = 1,nnodes_r

                irow = (iz-1)*nnodes_r + ir ! Current row

                r = (dble(ir)-0.5d0)*dr
                sr2 = diffusivity_r * dt / (2.0d0*r*dr)

                if (iz.eq.1) then                                 ! axial symmetry BC
                    a(i) = 1.0d0
                    row_indices(i) = irow
                    i = i + 1
                elseif (ir.eq.1) then                             ! radial symmetry BC
                    a(i) = 1.0d0
                    row_indices(i) = irow
                    i = i + 1
                elseif (iz.eq.nnodes_z) then                      ! constant BC
                    a(i) = 1.0d0
                    row_indices(i) = irow
                    i = i + 1
                elseif (ir.eq.nnodes_z) then                      ! constant BC
                    a(i) = 1.0d0
                    row_indices(i) = irow
                    i = i + 1                    
                elseif (icol.eq.irow-1) then
                    a(i) = beta*(sr2-sr1)                         ! u(i-1,j  ,k+1) coefficient
                    row_indices(i) = irow
                    i = i + 1
                elseif (icol.eq.irow) then
                    a(i) = 1.0d0 + 2.0d0*sz*beta + 2.0d0*sr1*beta ! u(i  ,j  ,k+1) coefficient
                    row_indices(i) = irow
                    i = i + 1
                elseif (icol.eq.irow+1) then
                    a(i) = -beta*(sr2+sr1)                        ! u(i-1,j  ,k+1) coefficient
                    row_indices(i) = irow
                    i = i + 1
                elseif (icol.eq.(iz-2)*nnodes_r+ir) then
                    a(i) = -sz*beta                               ! u(i  ,j-1,k+1) coefficient
                    row_indices(i) = irow
                    i = i + 1
                elseif (icol.eq.(iz-0)*nnodes_r+ir) then
                    a(i) = -sz*beta                               ! u(i  ,j+1,k+1) coefficient
                    row_indices(i) = irow
                    i = i + 1
                endif

            enddo
        enddo

    enddo
    column_pointers(n+1) = nnz+1



    ! Load b matrix
    iz = 1
    do ir = 1,nnodes_r
        irow = (iz-1)*nnodes_r + ir ! Current row
        b(irow) = u(ir,iz)
    enddo
    iz = nnodes_z
    do ir = 1,nnodes_r
        irow = (iz-1)*nnodes_r + ir ! Current row
        b(irow) = u(ir,iz)
    enddo
    ir = 1
    do iz = 1,nnodes_z
        irow = (iz-1)*nnodes_r + ir ! Current row
        b(irow) = u(ir,iz)
    enddo
    ir = nnodes_r
    do iz = 1,nnodes_z
        irow = (iz-1)*nnodes_r + ir ! Current row
        b(irow) = u(ir,iz)
    enddo
    do iz = 2,nnodes_z-1 ! Loop through all non-BC rows
        do ir = 2,nnodes_r-2
            irow = (iz-1)*nnodes_r + ir ! Current row
            if (ir.eq.1) then
                b(irow) = u(ir,iz)
            elseif (iz.eq.1) then
                b(irow) = u(ir,iz)
            elseif (ir.eq.nnodes_r.or.iz.eq.nnodes_z) then
                b(irow) = u(ir,iz)
            else
                b(irow) = (1.0d0 - 2.0d0*sz*beta - 2.0d0*sr1*beta)*u(ir,iz) &
                                     + sr1*(1.0d0-beta)*(u(ir-1,iz)+u(ir+1,iz)) &
                                     +  sz*(1.0d0-beta)*(u(ir,iz-1)+u(ir,iz+1)) &
                                     + sr2*(1.0d0-beta)*(u(ir+1,iz)-u(ir-1,iz))
            endif
        enddo
    enddo




    ! Per SuperLU instructions, solve the sparse matrix equation
    ! call c_fortran_dgssv(iopt,n,nnz,nrhs,values,row_indices,column_pointers,b,ldb,factors,info)
    call c_fortran_dgssv()


    return
    end subroutine

#endif


!--------------------------------------------------------------------------------------------------!




    subroutine tridiag(a,b,c,r,u,n)
    !----
    ! Numerical recipes tridiagonal matrix equation solver
    !----
    implicit none
    integer :: n
    double precision :: a(n), b(n), c(n), r(n), u(n)
    ! Solves for a vector u(1:n) of length n the tridiagonal linear set given by equation (2.4.1).
    ! a(1:n), b(1:n), c(1:n), and r(1:n) are input vectors and  are not modified.
    ! Parameter: NMAX is the maximum expected value of n.
    INTEGER ::  j
    double precision :: bet, gam(n) ! One vector of workspace, gam is needed.
    bet = b(1)
    u(1) = r(1)/bet
    do j = 2,n ! Decomposition and forward substitution.
        gam(j) = c(j-1)/bet
        bet = b(j)-a(j)*gam(j)
        if (abs(bet).lt.1.0d-8) then ! Algorithm fails
            write(0,*) 'tridiag failed'
            return
        endif
        u(j) = (r(j)-a(j)*u(j-1))/bet
    enddo
    do j = n-1,1,-1 ! Backsubstitution.
        u(j) = u(j)-gam(j+1)*u(j+1)
    enddo
    return
    end subroutine


end module
