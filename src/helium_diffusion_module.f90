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
! from Piotraschke (2012):
!
! Piotraschke, R. (2012). Thermal and Geologic Constraints on the Cretaceous-to-Neogene Tectonic
!     Development of the Klamath Mountains, Northern California. M.Sc. Thesis, Pennsylvania State
!     University.
!----

double precision :: mol_he4             ! moles of helium-4
double precision :: mol_th232           ! moles of thorium-232
double precision :: mol_u238            ! moles of uranium-238

!--------------------------------------------------------------------------------------------------!
contains
!--------------------------------------------------------------------------------------------------!

function apatite_he_age(mass_th232,mass_u238,vol_he4) result(age)
!----
! Calculate the (U-Th)/He age of an apatite grain
!
! Inputs:
!     - mass_th232: mass of thorium-232 (nanograms)
!     - mass_u238: mass of uranium-238 (nanograms)
!     - vol_he4: volume of helium-4 (nano-cubic-centimeters)
!
! Rachel Piotraschke's Matlab function heagecalc
!----

use radioactivity, only: atomic_mass_u238, &
                         atomic_mass_th232
use fmm, only: zeroin

implicit none

! Arguments
double precision :: mass_th232, mass_u238, vol_he4, age

! Calculate moles of helium-4, thorium-232, and uranium-238
mol_he4 = vol_he4/2.24132d13
mol_th232 = mass_th232/atomic_mass_th232/1.0d9
mol_u238 = mass_u238/atomic_mass_u238/1.0d9

! Calculate age by finding zeros of age function (apatite_he_age_fn)
age = zeroin(0.0d0,1000.0d0,apatite_he_age_fn,1.0d-10)

return
end function

!--------------------------------------------------------------------------------------------------!

function apatite_he_age_fn(t) result(f)
!----
! The zeros of this radioactive decay function are the apatite-He age
!----

use radioactivity, only: decay_u238, &
                         decay_u235, &
                         decay_th232

implicit none

! Arguments
double precision :: f
double precision, intent(in) :: t

f = 8.0d0*mol_u238*(exp(decay_u238*t)-1.0d0) &
        + 7.0d0*mol_u238/137.88d0*(exp(decay_u235*t)-1.0d0) &
        + 6.0d0*mol_th232*(exp(decay_th232*t)-1.0d0) &
        - mol_he4

return
end function

!--------------------------------------------------------------------------------------------------!

subroutine apatite_he_conc(temp,time,n,he_conc_surf,he_conc_init,sphere_radius)
!----
! Given a temperature-time history, calculate the expected concentration of helium in an apatite
! grain over time
!----

use radioactivity, only: atomic_mass_u238, atomic_mass_th232, decay_u238, decay_u235, decay_th232

implicit none

! Arguments
integer :: n                            ! length of time and temperature arrays
double precision :: temp(n)             ! temperature array (C)
double precision :: time(n)             ! time array (Ma)
double precision :: he_conc_surf        ! surface helium concentration
double precision :: he_conc_init        ! initial helium concentration
double precision :: sphere_radius       ! radius of the apatite grain (microns)

! Local variables
integer :: i, k
double precision :: dt_ma, dt
double precision :: dr_microns, dr
integer :: nnodes
integer :: ntimesteps
double precision, allocatable :: radius(:), volume(:)
double precision, allocatable :: he_conc(:), he_conc_sub(:), he_conc_sub_new(:), he_mol(:), he_total
double precision, allocatable :: diffusivity(:)
double precision :: diffusion_number
double precision :: mass_th232_ng, mass_u238_ng, mass_th232, mass_u238, mol_u238, mol_th232
double precision, allocatable :: mol_conc_u238(:), mol_conc_th232(:)
double precision :: volume_grain
double precision, allocatable :: production_rate(:)

double precision, parameter :: diffusivity0 = 5.0d-3   ! Farley (2000) (m^2/s)
double precision, parameter :: activation_energy = 32.9d0 ! Farley (2000) (kcal/mol)
double precision, parameter :: gas_constant = 1.9858775d-3 ! Gas constant (kcal/mol/K)

double precision :: beta

double precision, allocatable :: a(:,:)
double precision, allocatable :: b(:)
! double precision :: mol_u238, mol_th232


beta = 0.75d0
dt_ma = 0.0005d0                            ! Time step size (Ma)
dt = dt_ma*3.15d13                          ! Time step size (seconds)
dr_microns = 0.5d0                          ! Spatial step size (microns)
dr = dr_microns*1.0d-6                      ! Spatial step size (meters)
nnodes = int(sphere_radius/dr_microns)+2    ! Number of spatial nodes
ntimesteps = int(time(n)/dt_ma)             ! Number of timesteps

mass_th232_ng = 0.1d0                       ! Mass of thorium-232 (nanograms)
mass_u238_ng = 0.1d0                        ! Mass of uranium-238 (nanograms)

allocate(radius(nnodes))
allocate(volume(nnodes))
allocate(he_conc(nnodes))
allocate(he_conc_sub(nnodes))
allocate(he_conc_sub_new(nnodes))
allocate(he_mol(nnodes))
allocate(diffusivity(ntimesteps))
allocate(mol_conc_th232(nnodes))
allocate(mol_conc_u238(nnodes))
allocate(production_rate(nnodes-2))
allocate(a(nnodes-2,3))
allocate(b(nnodes-2))


! Radius at each spatial node (meters)
do i = 1,nnodes
    radius(i) = (dble(i)-1.5d0)*dr
enddo

! Volume of each spherical shell
volume(1) = 0.0d0
volume(2) = 4.0d0/3.0d0*3.14159265d0*radius(2)**3
do i = 3,nnodes
    volume(i) = 4.0d0/3.0d0*3.14159265d0*(radius(i)**3-radius(i-1)**3)
enddo


! Diffusivity history (function of temperature)
! NOTE: WILL NEED TO MAP INPUT TEMPERATURE TIMING TO NEW TIMING
diffusivity = 0.0d0
do i = 1,ntimesteps
    diffusivity(i) = diffusivity0*exp(-activation_energy/(gas_constant*(temp(i)+273.15d0)))
enddo


! Initial and boundary conditions
he_conc = he_conc_init
he_conc(nnodes) = he_conc_surf
do i = 1,nnodes
    he_conc_sub(i) = he_conc(i)*radius(i)
enddo
he_conc_sub(1) = -he_conc_sub(2) ! NOTE: SHOULD THIS BE POSITIVE?


! Helium production
! Rachel's note: this is constant in time but could vary in space if uranium and thorium are zoned
! Mass of uranium-238 and thorium-232 in the grain
mass_th232 = mass_th232_ng*1.0d-9
mass_u238 = mass_u238_ng*1.0d-9

! Moles of uranium-238 and thorium-232 in the grain
mol_th232 = mass_th232/atomic_mass_th232
mol_u238 = mass_u238/atomic_mass_u238

! Grain volume (m^3)
volume_grain = 4.0d0/3.0d0*3.14159265d0*(sphere_radius*1.0d-6)**3

! Molar concentration at each node (mol/m^3)
mol_conc_th232(2:nnodes-1) = mol_th232/volume_grain
mol_conc_u238(2:nnodes-1) = mol_u238/volume_grain

! Production rate at each node (mol/m^3/s)
production_rate(1) = 0.0d0
production_rate(nnodes) = 0.0d0
do i = 2,nnodes-1
    production_rate(i) = 8.0d0*mol_conc_u238(i)*decay_u238 + &
                         7.0d0*mol_conc_u238(i)/137.88d0*decay_u235 + &
                         6.0d0*mol_conc_th232(i)*decay_th232
enddo


! Step through time and calculate helium concentration
do k = 2,ntimesteps
    diffusion_number = diffusivity(k)*dt/dr**2

    if (temp(k).gt.90.0d0) then
        he_conc_sub_new = 0.0d0         ! All helium escapes above 90 C

    else
        ! Update the helium concentration at all nodes with newly produced radiogenic helium
        he_conc(nnodes) = he_conc_surf
        do i = 2,nnodes-1
            he_conc(i) = he_conc(i) + production_rate(i)*dt
        enddo

        ! Make the substitution
        do i = 1,nnodes
            he_conc_sub(i) = he_conc(i)*radius(i)
        enddo
        he_conc_sub(1) = -he_conc_sub(2)

        ! Load the tridiagonal diffusion matrix
        do i = 1,nnodes-2
            a(i,1) = -diffusion_number*beta
            a(i,2) = 2.0d0*diffusion_number*beta + 1.0d0
            a(i,3) = -diffusion_number*beta
        enddo
        a(1,2) = 1.0d0;
        a(1,3) = 0.0d0;
        a(nnodes-2,nnodes-3) = 0.0d0;
        a(nnodes-2,nnodes-2) = 1.0d0;

        ! Load the right-hand-side vector
        do i = 1,nnodes-2
            b(i) = diffusion_number*(1.0d0-beta)*he_conc_sub(i) + &
                   (1.0d0-(2.0d0*diffusion_number*(1.0d0-beta)))*he_conc_sub(i+1) + &
                   diffusion_number*(1.0d0-beta)*he_conc_sub(i+2)
        enddo

        ! Solve for the new helium concentration (the substituted function)
        call tridiag(a(:,1),a(:,2),a(:,3),b,he_conc_sub_new(2:nnodes-1),nnodes-2)
    endif

    ! Update the helium concentration
    do i = 2,nnodes-1
        he_conc(i) = he_conc_sub_new(i)/radius(i)
    enddo

    ! Moles of He in each shell
    do i = 2,nnodes-1
        he_mol(i) = volume(i)*he_conc(i)
    enddo

    ! Total He in the grain
    he_total = 0.0d0
    do i = 2,nnodes-1
        he_total = he_total + he_mol(i) ! moles
    enddo
    he_total = he_total*2.24132d13 ! nano-cubic-centimeters

enddo


return
end subroutine

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



end module !=======================================================================================!







program main
use apatite_he, only: apatite_he_age
implicit none
print *,apatite_he_age(1.0d0,1.0d0,0.05d0)
end
!
!
!
!
!
!
!
!
!
!
!
!
!
!
! subroutine apatite_he_2()
! implicit none
!
!
! end
