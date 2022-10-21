module helium_diffusion
!----
! Helium diffusion module
!
! This module contains procedures for studying the thermochronology of radiogenic helium
!
! Determine the amount of radiogenic helium and the corresponding thermochronologic age in
! grains that have undergone some thermal history by solving the production-diffusion equation
! in a spherical coordinate system (Carslaw and Jaeger, 1959):
!
!   dH            1       d           dH
!  ----  =  D * ----- * ---- [ r^2 * ---- ]  +  P
!   dt           r^2     dr           dr
!
! where
!     - H: Helium concentration
!     - P: Rate of radiogenic helium production
!     - D: Diffusivity of helium (thermally controlled)
!     - r: distance from the center of the grain
!
! The derivation of the finite difference approximation and the original Matlab program HeAge are
! from Piotraschke (2012).
!
! References:
!     Farley, K.A. (2000). Helium diffusion from apatite: General behavior as illustrated by
!         Durango fluorapatite. Journal of Geophysical Research, 105, B2, 2903-2914.
!     Piotraschke, R. (2012). Thermal and Geologic Constraints on the Cretaceous-to-Neogene Tectonic
!         Development of the Klamath Mountains, Northern California. M.Sc. Thesis, Pennsylvania State
!         University.
!----

! Physical constants
double precision, parameter :: univ_gas_constant = 1.9858775d-3             ! [kcal/mol/K]        R

! Atomic constants (IUGS 1977)
double precision, parameter :: atomic_mass_u238 = 238.0507826d0             !                     amU238
double precision, parameter :: atomic_mass_th232 = 232.0380553d0            !                     amTh232
double precision, parameter :: decay_u238 = 1.55125d-4                      ! [1/Ma]
double precision, parameter :: decay_u235 = 9.8485d-4                       ! [1/Ma]
double precision, parameter :: decay_th232 = 4.9475d-5                      ! [1/Ma]


! Diffusivity of helium in apatite (Farley, 2000: p. 2910)
! Note: this ranges from 0.0008 for c-perpendicular diffusion to 0.0130 for c-parallel diffusion
double precision, parameter :: he_diffusivity_apatite = 0.005d0             ! [m^2/s]             Do

! Helium diffusion activation energy in apatite (Farley, 2000: p. 2910)
double precision, parameter :: he_activation_energy_apatite = 32.9d0        ! [kcal/mol]          Ea


! Spherical finite difference parameters (generated in init_spherical_node_geometry)
integer :: nnodes_sphere
double precision, allocatable :: radius_shell(:)                            ! [m]
double precision, allocatable :: volume_shell(:)                            ! [m^3]
double precision :: volume_grain                                            ! [m^3]


!==================================================================================================!
contains
!==================================================================================================!

! ! function apatite_he_age(mass_th232,mass_u238,vol_he4) result(age)
! ! !----
! ! ! Calculate the (U-Th)/He age of an apatite grain
! ! !
! ! ! Inputs:
! ! !     - mass_th232: mass of thorium-232 (nanograms)
! ! !     - mass_u238: mass of uranium-238 (nanograms)
! ! !     - vol_he4: volume of helium-4 (nano-cubic-centimeters)
! ! !
! ! ! Rachel Piotraschke's Matlab function heagecalc
! ! !----

! use fmm, only: zeroin

! ! implicit none

! ! ! Arguments
! ! double precision :: mass_th232, mass_u238, vol_he4, age

! ! ! Calculate moles of helium-4, thorium-232, and uranium-238
! ! mol_he4 = vol_he4/2.24132d13
! ! mol_th232 = mass_th232/atomic_mass_th232/1.0d9
! ! mol_u238 = mass_u238/atomic_mass_u238/1.0d9

! ! ! Calculate age by finding zeros of age function (apatite_he_age_fn)
! ! age = zeroin(0.0d0,1000.0d0,apatite_he_age_fn,1.0d-10)

! ! return
! ! end function

! !--------------------------------------------------------------------------------------------------!

! ! function apatite_he_age_fn(t) result(f)
! ! !----
! ! ! The zeros of this radioactive decay function are the apatite-He age
! ! !----

! ! use radioactivity, only: decay_u238, &
! !                          decay_u235, &
! !                          decay_th232

! ! implicit none

! ! ! Arguments
! ! double precision :: f
! ! double precision, intent(in) :: t

! ! f = 8.0d0*mol_u238*(exp(decay_u238*t)-1.0d0) &
! !         + 7.0d0*mol_u238/137.88d0*(exp(decay_u235*t)-1.0d0) &
! !         + 6.0d0*mol_th232*(exp(decay_th232*t)-1.0d0) &
! !         - mol_he4

! ! return
! ! end function

!--------------------------------------------------------------------------------------------------!

subroutine calc_he_conc_apatite(temp_celsius, n, dt_ma, radius_microns)
!----
! Given a temperature-time history, calculate the expected concentration of helium in an apatite
! grain over time
!
! Inputs:
!   n:                number of time steps
!   temp_celsius:     n-point double precision array of temperatures (Celsius)
!   dt_ma:            time step size (Ma)
!   radius_microns:   radius of the (spherical) apatite grain (microns)
!----

implicit none

! Arguments
integer :: n
double precision :: temp_celsius(n)
double precision :: dt_ma
double precision :: radius_microns

! Local variables
double precision :: temp_kelvin(n)
double precision :: arg
integer :: i
double precision :: he_conc_surf        ! surface helium concentration
double precision :: he_conc_init        ! initial helium concentration
double precision :: dt_seconds
double precision :: dr_microns
double precision :: dr_meters
double precision, allocatable :: he_conc(:) 
double precision, allocatable :: he_conc_sub(:) 
double precision, allocatable :: he_conc_sub_new(:)
double precision, allocatable :: he_mol(:)
double precision, allocatable :: he_total(:)
double precision, allocatable :: diffusivity(:)
double precision :: diffusion_number
double precision, allocatable :: production_rate(:)

integer :: itime

double precision :: beta

double precision, allocatable :: a(:,:)
double precision, allocatable :: b(:)



! Convert temperature and time units
temp_kelvin = temp_celsius + 273.0d0
dt_seconds = dt_ma*1.0d6*365.25d0*24d0*60d0*60d0

beta = 0.85d0


!----
! Set up finite difference parameters (assuming spherical apatite grain)
!----
! Generate spherical nodes and geometric variables:
!     - nnodes_sphere
!     - radius_shell
!     - volume_shell
dr_microns = 0.5d0                ! Spatial step size (microns)
dr_meters = dr_microns*1.0d-6     ! Spatial step size (meters)
call init_spherical_node_geometry(radius_microns,dr_microns)

! Allocate helium concentration arrays
if (.not.allocated(he_conc)) then
    allocate(he_conc(nnodes_sphere))
endif
if (.not.allocated(he_conc_sub)) then
    allocate(he_conc_sub(nnodes_sphere))
endif
if (.not.allocated(he_conc_sub_new)) then
    allocate(he_conc_sub_new(nnodes_sphere))
endif
if (.not.allocated(he_mol)) then
    allocate(he_mol(nnodes_sphere))
endif

! Initialize helium concentration values and boundary conditions
he_conc_init = 0.0d0          ! Initial helium concentration
he_conc_surf = 0.0d0          ! Surface helium concentration
he_conc = he_conc_init        ! Set all nodes to have initial helium concentration
he_conc(nnodes_sphere) = he_conc_surf     ! Except surface...

! Set radiogenic helium production rate
if(.not.allocated(production_rate)) then
    allocate(production_rate(nnodes_sphere))
endif
call set_he_production_rate_apatite(production_rate,nnodes_sphere)

! Diffusivity history at each time step (Farley, 2000)
if (.not.allocated(diffusivity)) then
    allocate(diffusivity(n))
endif
diffusivity = 0.0d0
do i = 1,n
    arg = -he_activation_energy_apatite/(univ_gas_constant*temp_kelvin(i))
    diffusivity(i) = he_diffusivity_apatite * exp(arg)
enddo




!----
! Determine helium concentration throughout apatite grain over temperature history
!----
! Allocate and initialize array for total helium in grain (nano-cubic-centimeters)
allocate(he_total(n))
he_total = 0.0d0

! Allocate arrays for solving helium concentration at each time step
allocate(a(nnodes_sphere-2,3))
allocate(b(nnodes_sphere-2))

! Step through time and calculate helium concentration
do itime = 2,n

    ! Calculate diffusion number
    diffusion_number = diffusivity(itime)*dt_seconds/dr_meters**2

!     if (temp(k).gt.90.0d0) then
!         he_conc_sub_new = 0.0d0         ! All helium escapes above 90 C

!     else
        ! Update the helium concentration at all nodes with newly produced radiogenic helium
        he_conc(nnodes_sphere) = he_conc_surf
        do i = 2,nnodes_sphere-1
            he_conc(i) = he_conc(i) + production_rate(i)*dt_seconds
        enddo

        ! Make the substitution for spherical coordinates
        do i = 1,nnodes_sphere
            he_conc_sub(i) = he_conc(i)*radius_shell(i)
        enddo
        he_conc_sub(1) = -he_conc_sub(2)

        ! Load the tridiagonal diffusion matrix with the current diffusion number
        do i = 1,nnodes_sphere-2
            a(i,1) = -diffusion_number*beta
            a(i,2) = 2.0d0*diffusion_number*beta + 1.0d0
            a(i,3) = -diffusion_number*beta
        enddo
        a(1,2) = 1.0d0;
        a(1,3) = 0.0d0;
        a(nnodes_sphere-2,1) = 0.0d0;
        a(nnodes_sphere-2,2) = 1.0d0;

        ! Load the right-hand-side vector with diffusion number and helium concentration
        do i = 1,nnodes_sphere-2
            b(i) = diffusion_number*(1.0d0-beta)*he_conc_sub(i) + &
                   (1.0d0-(2.0d0*diffusion_number*(1.0d0-beta)))*he_conc_sub(i+1) + &
                   diffusion_number*(1.0d0-beta)*he_conc_sub(i+2)
        enddo

        ! Solve for the new helium concentration (the substituted function)
        he_conc_sub_new = he_conc_sub
        call tridiag(a(:,1),a(:,2),a(:,3),b,he_conc_sub_new(2:nnodes_sphere-1),nnodes_sphere-2)
!     endif

    ! Update the helium concentration
    he_conc = he_conc_sub_new/radius_microns ! SHOULD THIS BE DIVIDED BY R AT EACH NODE???

    ! Moles of He in each shell
    do i = 2,nnodes_sphere-1
        he_mol(i) = volume_shell(i)*he_conc(i)
    enddo

    ! Total He in the grain
    he_total(itime) = 0.0d0
    do i = 2,nnodes_sphere-1
        he_total(itime) = he_total(itime) + he_mol(i) ! moles
    enddo
    he_total(itime) = he_total(itime)*2.24132d13 ! nano-cubic-centimeters

enddo


return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine init_spherical_node_geometry(dr_microns,radius_microns)
!----
! Create spherical node geometry for finite difference procedure
!
!
!         dr_microns                                \
!          -------                                   |
!          *     *     *     *     *     *     *     * radius_microns
!                                                    |
!                                                   /
!
! Inputs:
!   dr_microns:       initial track length (microns)
!   radius_microns:   number of time steps
!
! Outputs (helium_diffusion module variables)
!   nnodes_sphere:    number of finite difference nodes
!   radius_shell:     distance from center of sphere to each node [m]
!   volume_shell:     volume of shell corresponding to each node  [m^3]
!   volume_sphere:    volume of entire sphere                     [m^3]
!----

implicit none


! Arguments
double precision :: dr_microns
double precision :: radius_microns

! Local variables
integer :: i


! Number of spatial nodes + boundary condition nodes
nnodes_sphere = int(radius_microns/dr_microns)+2    ! Number of spatial nodes

! Allocate memory to shell radius/volume arrays
if (.not.allocated(radius_shell)) then
    allocate(radius_shell(nnodes_sphere))
endif
if (.not.allocated(volume_shell)) then
    allocate(volume_shell(nnodes_sphere))
endif

! Distance from center of grain at each spatial node (meters)
do i = 1,nnodes_sphere
    radius_shell(i) = (dble(i)-1.5d0)*dr_microns*1.0d-6
enddo

! Volume of each spherical shell
volume_shell(1) = 0.0d0
volume_shell(2) = 4.0d0/3.0d0*3.14159265d0*radius_shell(2)**3
do i = 3,nnodes_sphere
    volume_shell(i) = 4.0d0/3.0d0*3.14159265d0*(radius_shell(i)**3-radius_shell(i-1)**3)
enddo

! Sphere volume (m^3)
volume_grain = 4.0d0/3.0d0*3.14159265d0*(radius_microns*1.0d-6)**3

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine set_he_production_rate_apatite(p,n)
!----
! Set radiogenic helium production rate at all nodes of a spherical apatite grain by assuming an
! initial uranium and thorium concentration
!
! Rachel's note: this is constant in time but could vary in space if uranium and thorium are zoned
!----

implicit none

! Arguments
integer :: n
double precision :: p(n)

! Local variables
double precision :: mass_u238_ng
double precision :: mass_th232_ng
double precision :: mass_u238_grams
double precision :: mass_th232_grams
double precision :: mol_u238
double precision :: mol_th232
double precision :: mol_conc_u238(n)
double precision :: mol_conc_th232(n)
integer :: i


! Mass of uranium-238 and thorium-232 in the grain
mass_u238_ng = 0.1d0                        ! Initial mass of uranium-238 (nanograms)
mass_th232_ng = 0.1d0                       ! Initial mass of thorium-232 (nanograms)
mass_u238_grams = mass_u238_ng*1.0d-9       ! Initial mass of uranium-238 (grams)
mass_th232_grams = mass_th232_ng*1.0d-9     ! Initial mass of thorium-232 (nanograms)

! Moles of uranium-238 and thorium-232 in the grain
mol_u238 = mass_u238_grams/atomic_mass_u238
mol_th232 = mass_th232_grams/atomic_mass_th232

! Molar concentration at each node (mol/m^3)
mol_conc_u238 = mol_u238/volume_grain
mol_conc_th232 = mol_th232/volume_grain

! Production rate at each node (mol/m^3/s)
p(1) = 0.0d0
p(n) = 0.0d0
do i = 2,n-1
    p(i) = 8.0d0*mol_conc_u238(i)*decay_u238 + &
           7.0d0*mol_conc_u238(i)/137.88d0*decay_u235 + &
           6.0d0*mol_conc_th232(i)*decay_th232
enddo

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






