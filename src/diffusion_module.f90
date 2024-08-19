!--------------------------------------------------------------------------------------------------!
! Module Diffusion                                                                                 !
!                                                                                                  !
! Authors:                                                                                         !
!     - Matt Herman (current Modern Fortran version, i.e., what you are looking at right now!)     !
!     - Rachel Piotraschke (original version; PSU MS thesis)                                       !
!     - Kevin Furlong (original version; supervisor)                                               !
!                                                                                                  !
! This module contains variables and subroutines numerically solving the production-diffusion      !
! differential equation. The derivation of the finite difference approximation and the original    !
! Matlab program HeAge are from Piotraschke (2012) (MS thesis).                                    !
!                                                                                                  !
! The (production-)diffusion equation in a Cartesian coordinate system can be written as (Carslaw  !
! and Jaeger, 1959):                                                                               !
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
!     u: function that will undergo diffusion (temperature, concentration, topography, etc.)       !
!     t: time                                                                                      !
!     x: position                                                                                  !
!     r: radius                                                                                    !
!     D: diffusivity                                                                               !
!     P: rate of production                                                                        !
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