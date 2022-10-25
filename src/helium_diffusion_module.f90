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
double precision, parameter :: decay_u238 = 1.55125d-4/3.15d13              ! [1/s]
double precision, parameter :: decay_u235 = 9.8485d-4/3.15d13               ! [1/s]
double precision, parameter :: decay_th232 = 4.9475d-5/3.15d13              ! [1/s]


! Diffusivity of helium in apatite (Farley, 2000: p. 2910)
! Note: this ranges from 0.0008 for c-perpendicular diffusion to 0.0130 for c-parallel diffusion
double precision, parameter :: he_diffusivity_apatite = 0.005d0             ! [m^2/s]             Do

! Helium diffusion activation energy in apatite (Farley, 2000: p. 2910)
double precision, parameter :: he_activation_energy_apatite = 32.9d0        ! [kcal/mol]          Ea


! Spherical finite difference parameters (generated in init_spherical_node_geometry)
integer :: nnodes_sphere
double precision, allocatable :: radius_shell(:)                            ! [m]
double precision, allocatable :: volume_shell(:)                            ! [m^3]
double precision :: volume_sphere                                           ! [m^3]


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



subroutine calc_he_conc_apatite(temp_celsius, n, dt_ma, radius_microns, beta)
!----
! Given a temperature-time history, calculate the concentration of helium in an a spherical apatite
! grain over time using a finite difference approximation for the diffusion equation
!
! Inputs:
!   n:                number of time steps
!   temp_celsius:     n-point double precision array of temperatures (Celsius)
!   dt_ma:            time step size for the temperature history (Ma)
!   radius_microns:   radius of the (spherical) apatite grain (microns)
!   beta:             implicitness weight in finite difference solution (0-1)
!----

implicit none

! Arguments
integer :: n
double precision :: temp_celsius(n)
double precision :: dt_ma
double precision :: radius_microns
double precision :: beta

! Local variables
double precision, allocatable :: temp_kelvin(:)
double precision :: temp_max_kelvin
double precision, allocatable :: diffusivity(:)
double precision :: diffusion_number
double precision :: diffusivity_max
double precision :: dt_seconds
double precision :: dt_seconds_resamp
integer :: nresamp
double precision :: dr_microns
double precision :: dr_meters
double precision :: he_conc_surf
double precision :: he_conc_init
double precision, allocatable :: he_conc(:) 
double precision, allocatable :: he_conc_sub(:) 
double precision, allocatable :: he_conc_sub_new(:)
double precision, allocatable :: he_mol(:)
double precision, allocatable :: he_total(:)
double precision, allocatable :: production_rate(:)
double precision :: arg
integer :: i, j
integer :: itime
double precision, allocatable :: a(:,:)
double precision, allocatable :: b(:,:)



!----
! Set up finite difference spatial grid (assuming spherical apatite grain)
!----
write(*,*) 'calc_he_conc_apatite: setting up finite difference spatial grid' 

! Generate spherical nodes and geometric variables (from the helium_diffusion module)
!     - nnodes_sphere           spatial nodes + 2 BC nodes (at sphere center and beyond sphere edge)
!     - radius_shell [m]        distance of each node from center of sphere
!     - volume_shell [m^3]      volume of spherical shell for each node
!     - volume_sphere [m^3]     volume of entire sphere
nnodes_sphere = 102
dr_microns = radius_microns/dble(nnodes_sphere-2) ! Spatial step size (microns)
dr_meters = dr_microns*1.0d-6                     ! Spatial step size (meters)
call init_spherical_node_geometry(radius_microns,dr_microns)
write(*,*) 'radius_microns=',radius_microns
write(*,*) 'dr_microns=    ',dr_microns
write(*,*) 'nnodes_sphere= ',nnodes_sphere



!----
! Set up finite difference time stepping
!----
write(*,*) 'calc_he_conc_apatite: resampling temperature history onto new time grid' 

! Calculate time step size in seconds for input thermal history
dt_seconds = dt_ma*1.0d6*365.25d0*24.0d0*60.0d0*60.0d0

! Determine new time step size based on maximum diffusivity
temp_max_kelvin = maxval(temp_celsius) + 273.0d0
arg = -he_activation_energy_apatite/(univ_gas_constant*temp_max_kelvin)
diffusivity_max = he_diffusivity_apatite * exp(arg)
dt_seconds_resamp = 0.5d0*dr_meters**2/diffusivity_max ! Stability criterion for fully explicit problems

! Check resampling time step size
if (dt_seconds_resamp.lt.0.01d0*dt_seconds) then
    ! Implicit solution should be stable, so limit shrinking of resampled time step
    dt_seconds_resamp = 0.01d0*dt_seconds
elseif (dt_seconds_resamp.gt.dt_seconds) then
    ! Keep thermal history time step if smaller than resampled time step
    dt_seconds_resamp = dt_seconds
endif
write(*,*) 'dt_ma=',dt_ma
write(*,*) 'dt_seconds=',dt_seconds
write(*,*) 'dt_ma_resamp=',dt_seconds_resamp

! Calculate number of resampled time steps
nresamp = int(dble(n)*dt_seconds/dt_seconds_resamp)
write(*,*) 'n=',n
write(*,*) 'nresamp=',nresamp

! Resample thermal history onto new time grid
if (.not.allocated(temp_kelvin)) then
    allocate(temp_kelvin(nresamp))
endif
temp_kelvin = 0.0d0
do i = 1,nresamp
    j = int(dt_seconds_resamp*dble(i)/dt_seconds)
    if (j.lt.1) then
        j = 1
    endif
    temp_kelvin(i) = temp_celsius(j) + 273.0d0
enddo


!----
! Calculate diffusivity of helium in apatite over (resampled) temperature history
!----
write(*,*) 'calc_he_conc_apatite: calculating helium diffusivity over temperature history'

! Diffusivity at each time step (Farley, 2000) [m^2/s]
if (.not.allocated(diffusivity)) then
    allocate(diffusivity(nresamp))
endif
diffusivity = 0.0d0
do i = 1,nresamp
    arg = -he_activation_energy_apatite/(univ_gas_constant*temp_kelvin(i))
    diffusivity(i) = he_diffusivity_apatite * exp(arg)
enddo
write(*,*) 'min(diffusivity)=',minval(diffusivity)
write(*,*) 'max(diffusivity)=',maxval(diffusivity)



!----
! Calculate helium production rate
!----
write(*,*) 'calc_he_conc_apatite: calculating helium production rate'

! Set radiogenic helium production rate (mol/m^3/s)
if(.not.allocated(production_rate)) then
    allocate(production_rate(nnodes_sphere))
endif
call set_he_production_rate_apatite(production_rate,nnodes_sphere)



!----
! Initialize helium concentration
!----
write(*,*) 'calc_he_conc_apatite: initializing helium concentration'

! Allocate helium arrays
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
if (.not.allocated(he_total)) then
    allocate(he_total(nresamp))
endif

! Initialize helium concentration values and boundary conditions
he_conc_init = 1.0d0          ! Initial helium concentration
he_conc_surf = 0.0d0          ! Surface helium concentration
he_conc = he_conc_init        ! Set all nodes to have initial helium concentration
he_conc(nnodes_sphere) = he_conc_surf     ! Except surface boundary condition node...
write(*,*) 'he_conc_init=',he_conc_init
write(*,*) 'he_conc_surf=',he_conc_surf





!----
! Determine helium concentration in the apatite grain over temperature history
!----
open(unit=63,file='junk.out',status='unknown')

! Initialize total helium to be zero
he_total = 0.0d0


! Allocate arrays for determining helium concentration at each node at each time step
if (.not.allocated(a)) then
    allocate(a(nnodes_sphere,3))
endif
if (.not.allocated(b)) then
    allocate(b(nnodes_sphere,1))
endif


! Calculate helium concentration at each time step by solving finite difference
! approximation to the diffusion equation
do itime = 2,nresamp

    ! write(*,*) 'calc_he_conc_apatite: working on itime=',itime,' of',nresamp

    ! Calculate diffusion number
    diffusion_number = diffusivity(itime)*dt_seconds/dr_meters**2
    ! write(*,*) 'diffusion_number=',diffusion_number

!     if (temp(k).gt.90.0d0) then
!         he_conc_sub_new = 0.0d0         ! All helium escapes above 90 C

!     else
        ! ! Update the helium concentration at all nodes with newly produced radiogenic helium
        ! ! write(*,*) 'calc_he_conc_apatite: producing radiogenic helium'
        ! do i = 1,nnodes_sphere
        !    he_conc(i) = he_conc(i) + production_rate(i)*dt_seconds
        ! enddo

        ! Make the substitution for spherical coordinates
        ! write(*,*) 'calc_he_conc_apatite: converting helium concentration to spherical coordinates'
        do i = 1,nnodes_sphere
            he_conc_sub(i) = he_conc(i)*radius_shell(i)
        enddo
        he_conc_sub(1) = -he_conc_sub(2)

        ! Load the tridiagonal (implicit) diffusion matrix with the current diffusion number
        ! write(*,*) 'calc_he_conc_apatite: loading tridiagonal matrix'
        do i = 1,nnodes_sphere
            a(i,1) = -diffusion_number*beta
            a(i,2) = 2.0d0*diffusion_number*beta + 1.0d0
            a(i,3) = -diffusion_number*beta
        enddo
        a(1,1) = 0.0d0
        a(1,2) = 1.0d0
        a(1,3) = 0.0d0
        a(nnodes_sphere,1) = 0.0d0
        a(nnodes_sphere,2) = 1.0d0
        a(nnodes_sphere,3) = 0.0d0

        ! Load the right-hand-side vector with diffusion number and helium concentration
        ! write(*,*) 'calc_he_conc_apatite: loading right-hand-side vector'
        do i = 2,nnodes_sphere-1
            b(i,1) = diffusion_number*(1.0d0-beta)*he_conc_sub(i-1) + &
                   (1.0d0-(2.0d0*diffusion_number*(1.0d0-beta)))*he_conc_sub(i) + &
                   diffusion_number*(1.0d0-beta)*he_conc_sub(i+1)
        enddo
        b(1,1) = he_conc_sub(1)
        b(nnodes_sphere,1) = he_conc_sub(nnodes_sphere)

        ! Solve for the new helium concentration (the substituted function)
        ! write(*,*) 'calc_he_conc_apatite: solving equation for updated helium concentration'
        ! call tridiag(a(:,1),a(:,2),a(:,3),b(:,1),he_conc_sub_new(2:nnodes_sphere-1),nnodes_sphere-2)
        call dgtsv(nnodes_sphere, &
                   1, &
                   a(2:nnodes_sphere,1), &
                   a(:,2), &
                   a(1:nnodes_sphere-1,3), &
                   b, &
                   nnodes_sphere, &
                   i)
!     endif

    ! Update the helium concentration
    ! write(*,*) 'calc_he_conc_apatite: calculating helium concentration at each node'
    do i = 1,nnodes_sphere
        ! he_conc(i) = he_conc_sub_new(i)/radius_shell(i)
        he_conc(i) = b(i,1)/radius_shell(i)
        ! print *,he_conc(i),he_conc_sub_new(i),radius_shell(i)
    enddo
    he_conc(1) = he_conc(2)
    he_conc(nnodes_sphere) = he_conc_surf

    ! if (mod(itime,1000).eq.32) then
    if (itime.eq.2.or.itime.eq.1032) then
        write(63,'(A)') '>'
        do i = 1,nnodes_sphere
            write(63,*) radius_shell(i),he_conc(i)
        enddo
    endif

    ! ! Moles of He in each shell
    ! ! write(*,*) 'calc_he_conc_apatite: calculating moles of helium in each shell'
    ! do i = 2,nnodes_sphere-1
    !     he_mol(i) = volume_shell(i)*he_conc(i)
    ! enddo

    ! ! Total He in the grain
    ! ! write(*,*) 'calc_he_conc_apatite: calculating total helium in the grain'
    ! he_total(itime) = 0.0d0
    ! do i = 2,nnodes_sphere-1
    !     he_total(itime) = he_total(itime) + he_mol(i) ! moles
    ! enddo
    ! he_total(itime) = he_total(itime)*2.24132d13 ! nano-cubic-centimeters

enddo


return
end subroutine



!--------------------------------------------------------------------------------------------------!



subroutine init_spherical_node_geometry(radius_microns,dr_microns)
!----
! Create spherical node geometry for finite difference procedure
!
!   Center of sphere
!          |
!          |       dr_microns                          \ radius_microns
!          V        -------                             |
!       o  X  *     *     *     *     *     *     *     *     o
! BC node                                               |     BC node
!                                                      /
!
! Inputs (arguments):
!   radius_microns:   radius of sphere (microns)
!   dr_microns:       spatial step size (microns)
!
! Inputs (helium_diffusion module variables):
!   nnodes_sphere:    number of spatial nodes, including 2 BC nodes
!
! Outputs (helium_diffusion module variables)
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

! Volume of each spherical shell (m^3)
volume_shell(1) = 0.0d0
volume_shell(2) = 4.0d0/3.0d0*3.14159265d0*radius_shell(2)**3
do i = 3,nnodes_sphere
    volume_shell(i) = 4.0d0/3.0d0*3.14159265d0*(radius_shell(i)**3-radius_shell(i-1)**3)
enddo

! Sphere volume (m^3)
volume_sphere = 4.0d0/3.0d0*3.14159265d0*(radius_microns*1.0d-6)**3

return
end subroutine



!--------------------------------------------------------------------------------------------------!



subroutine set_he_production_rate_apatite(p,n)
!----
! Set radiogenic helium production rate at all nodes of a spherical apatite grain by assuming an
! initial uranium and thorium concentration
!
! Inputs (arguments):
!   n:   number of nodes
!
! Outputs (arguments)
!   p:   production rate at each node [mol/m^3/s] (NOTE: ARRAY MUST BE ALLOCATED ON INPUT!)
!
! Rachel's notes (Appendix A of her thesis):
!   - Distribution of parent nuclides within the grain is assumed to be homogeneous. Consequently,
!     any redistribution of helium within the grain due to α-ejection is assumed to cancel out
!     (Ketcham, 2005), and effects of a heterogeneous distribution of parent nuclides on the
!     alpha-ejection correction (Farley and Stockli, 2002) are neglected. The probability that the
!     actual distribution is homogeneous is unknown, although recent work by Farley et al. (2011)
!     suggests that U and Th enrichment at the grain rim is more likely in at least some samples;
!     however, since no measurements of U and Th concentration at a sub-grain scale were made as
!     part of this study, no better assumption can be made.
!   - Helium production is assumed to be constant in time, because variation in abundance of parent
!     nuclides is negligible over timescales of less than ~300 Ma (Wolf et al., 1998).
!   - This is constant in time but could vary in space if uranium and thorium are zoned
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
mass_th232_grams = mass_th232_ng*1.0d-9     ! Initial mass of thorium-232 (grams)

! Moles of uranium-238 and thorium-232 in the grain
mol_u238 = mass_u238_grams/atomic_mass_u238
mol_th232 = mass_th232_grams/atomic_mass_th232

! Molar concentration at each node (mol/m^3)
mol_conc_u238 = mol_u238/volume_sphere
mol_conc_th232 = mol_th232/volume_sphere

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






