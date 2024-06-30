module apatite_helium
!----
! This module contains procedures for studying the thermochronology of radiogenic helium in apatite
!
! The derivation of the finite difference approximation and the original Matlab program HeAge are
! from Piotraschke (2012). The diffusion routines are in diffusion_module.f90. The routines for
! analyzing (U-Th)/He ages are in radiogenic_helium_module.f90.
!
! References:
! Farley, K.A. (2000). Helium diffusion from apatite: General behavior as illustrated by Durango
!     fluorapatite. Journal of Geophysical Research, 105, B2, 2903-2914.
! Piotraschke, R. (2012). Thermal and Geologic Constraints on the Cretaceous-to-Neogene Tectonic
!     Development of the Klamath Mountains, Northern California. M.Sc. Thesis, Pennsylvania State
!     University.
!----

! Physical constants
double precision, parameter :: univ_gas_constant = 1.9858775d-3             ! [kcal/mol/K]        R

! Diffusivity of helium in apatite (Farley, 2000: p. 2910)
! Note: this ranges from 0.0008 for c-perpendicular diffusion to 0.0130 for c-parallel diffusion
double precision, parameter :: he_diffusivity_apatite = 0.005d0             ! [m^2/s]

! Helium diffusion activation energy in apatite (Farley, 2000: p. 2910)
double precision, parameter :: he_activation_energy_apatite = 32.9d0        ! [kcal/mol]

! Helium parameters (from calc_he_conc_apatite)
double precision, allocatable :: he_conc(:)                                 ! [mol/m^3/s]
double precision, allocatable :: he_total(:)                                ! [nano cm^3]


double precision, parameter :: ma2s = 1.0d6*365d0*24d0*60d0*60d0



!==================================================================================================!
contains
!==================================================================================================!



subroutine calc_apatite_he_age(temp_celsius, &
                               dep_km, &
                               n, &
                               dt_ma, &
                               dt_var, &
                               radius_microns, &
                               nnodes, &
                               beta, &
                               tau_max, &
                               age)
!----
! Given a temperature-time history, calculate the concentration of helium in an a spherical apatite
! grain over time using a finite difference approximation for the diffusion equation
!
! Inputs:
!   n:                number of time steps
!   temp_celsius:     n-point double precision array of temperatures (Celsius)
!   dep_km:           n-point double precision array of depths, positive up (km)
!   dt_ma:            time step size for the temperature history (Ma)
!   radius_microns:   radius of the (spherical) apatite grain (microns)
!   nnodes:           number of spatial nodes + 2 boundary condition nodes
!   beta:             implicitness weight in finite difference solution (0-1)
!   tau_max:          maximum dimensionless time to retain He and calculate diffusion (over 1 Ma)
!
! Output:
!   age:              age of grain (Ma)
!----

use diffusion, only: nnodes_sphere, &
                     radius_shell, &
                     volume_shell, &
                     volume_sphere, &
                     init_spherical_node_geometry, &
                     diffusion_step_dirichlet

use radiogenic_helium, only: atomic_mass_th232, &
                             atomic_mass_u238, &
                             decay_th232, &
                             decay_u235, &
                             decay_u238, &
                             calc_he_production_rate, &
                             calc_u_th_he_age

implicit none

! Arguments
integer :: n                                                ! number of input timesteps
double precision :: temp_celsius(n)                         ! temperature history [C]
double precision :: dep_km(n)                               ! depth history [km]
double precision :: dt_ma                                   ! input timestep size [Ma]
double precision :: dt_var                                  ! >0: dt reduction; <0: timestep [Ma]
double precision :: radius_microns                          ! grain radius [um]
integer :: nnodes                                           ! number of spatial nodes
double precision :: beta                                    ! implicitness parameter
double precision :: tau_max                                 ! Max dimensionless time to keep He
double precision :: age                                     ! (U-Th)/He age [Ma]

! Local variables
double precision :: radius_meters                           ! grain radius [m]
double precision :: dr_microns                              ! node spacing [microns]
double precision :: dr_meters                               ! node spacing [m]
double precision :: dt_seconds                              ! input time spacing [s]
double precision :: dt_seconds_resamp                       ! resampled time spacing [s]
integer :: nresamp                                          ! resampled number of time steps
double precision, allocatable :: temp_kelvin(:)             ! resampled temperature history [K]
double precision :: temp_max_kelvin                         ! maximum temperature [K]
double precision, allocatable :: diffusivity(:)             ! resampled diffusivity history [m^2/s]
double precision :: diffusivity_max                         ! maximum diffusivity [m^2/s]
double precision :: mass_th232                              ! Mass of thorium-232 [g]
double precision :: mass_u238                               ! Mass of uranium-238 [g]
double precision :: mol_th232                               ! Moles of thorium-232
double precision :: mol_u235                                ! Moles of uranium-235
double precision :: mol_u238                                ! Moles of uranium-238
double precision, allocatable :: production_rate(:)         ! radiogenic He production [mol/m^3/s]
double precision :: he_conc_surf                            ! He surface BC [mol/m^3/s]
double precision :: he_conc_init                            ! He initial BC [mol/m^3/s]
double precision, allocatable :: he_conc_sub(:)             ! radius x He concentration
double precision, allocatable :: he_mol(:)                  ! moles of He at each node [mol]
double precision :: diffusion_number                        ! diffusivity x dt / dr^2
double precision :: arg                                     ! dummy variable for exp arguments
double precision :: tau                                     ! Non-dimensional time
logical :: runDiffusion                                     ! Diffusion switch
integer :: istart                                           ! starting index for diffusion
integer :: i, j, itime                                      ! loop indices



!----
! Set up finite difference spatial grid (assuming spherical apatite grain)
!----
! write(*,*) 'calc_apatite_he_age: setting up finite difference spatial grid' 

! Generate spherical nodes and geometric variables (diffusion module):
!     - nnodes_sphere           spatial nodes + 2 BC nodes (at sphere center and beyond sphere edge)
!     - radius_shell [m]        distance of each node from center of sphere [m]
!     - volume_shell [m^3]      volume of spherical shell for each node [m^3]
!     - volume_sphere [m^3]     volume of entire sphere [m^3]
radius_meters = radius_microns/1.0d6             ! Radius of grain (meters)
dr_microns = radius_microns/dble(nnodes-2)        ! Spatial step size (microns)
dr_meters = dr_microns*1.0d-6                     ! Spatial step size (meters)
call init_spherical_node_geometry(radius_meters,dr_meters,nnodes)
! write(*,*) '    radius_meters=',radius_meters
! write(*,*) '    dr_meters=    ',dr_microns
! write(*,*) '    nnodes_sphere=',nnodes_sphere
! do i=1,nnodes_sphere
!     write(*,*) i,radius_shell(i), volume_shell(i)
! enddo
! write(*,*) '    volume_sphere=',volume_sphere



!----
! Set up finite difference time stepping
!----
! write(*,*) 'calc_apatite_he_age: resampling temperature history onto new time grid' 

! Calculate time step size in seconds for input thermal history
dt_seconds = dt_ma*ma2s

if (dt_var.gt.0.0d0) then
    ! Determine new time step size based on maximum diffusivity
    temp_max_kelvin = maxval(temp_celsius) + 273.0d0
    arg = -he_activation_energy_apatite/(univ_gas_constant*temp_max_kelvin)
    diffusivity_max = he_diffusivity_apatite * exp(arg)
    dt_seconds_resamp = 0.5d0*dr_meters**2/diffusivity_max ! Stability for fully explicit problems
    ! write(*,*) '    dt_ma_resamp=',dt_seconds_resamp

    ! Check resampling time step size
    if (dt_seconds_resamp.lt.dt_var*dt_seconds) then
        ! Implicit solution should be stable (IT IS NOT!), so limit shrinking of resampled time step
        dt_seconds_resamp = dt_var*dt_seconds
    elseif (dt_seconds_resamp.gt.dt_seconds) then
        ! Keep thermal history time step if smaller than resampled time step
        dt_seconds_resamp = dt_seconds
    endif
else
    dt_seconds_resamp = -dt_var*ma2s
endif
! write(*,*) '    dt_ma=       ',dt_ma
! write(*,*) '    dt_seconds=  ',dt_seconds
! write(*,*) '    dt_ma_resamp=',dt_seconds_resamp

! Calculate number of resampled time steps
nresamp = int(dble(n)*dt_seconds/dt_seconds_resamp)
! write(*,*) '    n=           ',n
! write(*,*) '    nresamp=     ',nresamp

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

! Determine starting time for diffusion calculation (rock has to exist!)
do i = 1,n
    if (dep_km(i).le.0.0d0) then ! Rock is born when it hits depth of 0
        istart = nint(dble(i)*dt_seconds/dt_seconds_resamp)
        if (istart.lt.1) then
            istart = 1
        endif
        exit
    endif
enddo


!----
! Calculate diffusivity of helium in apatite over (resampled) temperature history
!----
! write(*,*) 'calc_apatite_he_age: calculating helium diffusivity over temperature history'

! Diffusivity at each time step (Farley, 2000) [m^2/s]
if (.not.allocated(diffusivity)) then
    allocate(diffusivity(nresamp))
endif
diffusivity = 0.0d0
do i = 1,nresamp
    arg = -he_activation_energy_apatite/(univ_gas_constant*temp_kelvin(i))
    diffusivity(i) = he_diffusivity_apatite * exp(arg)
enddo
! write(*,*) '    min(diffusivity)=     ',minval(diffusivity)
! write(*,*) '    max(diffusivity)=     ',maxval(diffusivity)
! write(*,*) '    min(diffusion_number)=',minval(diffusivity)*dt_seconds_resamp/dr_meters**2
! write(*,*) '    max(diffusion_number)=',maxval(diffusivity)*dt_seconds_resamp/dr_meters**2


!----
! Calculate helium production rate
!----
! write(*,*) 'calc_apatite_he_age: calculating helium production rate'

! Calculate radiogenic helium production rate (mol/m^3/s), assuming uniform distribution of
! radioactive uranium and thorium
mass_th232 = 0.1d-9                       ! Initial mass of thorium-232 (grams)
mass_u238 = 0.1d-9                        ! Initial mass of uranium-238 (grams)
mol_th232 = mass_th232/atomic_mass_th232  ! Initial moles of thorium-232
mol_u238 = mass_u238/atomic_mass_u238     ! Initial moles of uranium-238
mol_u235 = mol_u238/137.88d0              ! Initial moles of uranium-235
if(.not.allocated(production_rate)) then
    allocate(production_rate(nnodes_sphere))
endif



!----
! Initialize helium concentration arrays
!----
! write(*,*) 'calc_apatite_he_age: initializing helium concentration'

! Allocate helium arrays
if (.not.allocated(he_conc)) then
    allocate(he_conc(nnodes_sphere))
endif
if (.not.allocated(he_conc_sub)) then
    allocate(he_conc_sub(nnodes_sphere))
endif
if (.not.allocated(he_mol)) then
    allocate(he_mol(nnodes_sphere))
endif
if (.not.allocated(he_total)) then
    allocate(he_total(nresamp))
endif

! Initialize helium concentration values and boundary conditions
he_conc_init = 0.0d0                      ! Initial helium concentration
he_conc_surf = 0.0d0                      ! Surface helium concentration
he_conc = he_conc_init                    ! Set all nodes to have initial helium concentration
he_conc(nnodes_sphere) = he_conc_surf     ! Except surface boundary condition node...
! write(*,*) '    he_conc_init=',he_conc_init
! write(*,*) '    he_conc_surf=',he_conc_surf





!----
! Determine helium concentration in the apatite grain over temperature history
!----
! write(*,*) 'calc_apatite_he_age: determining helium concentration over time'
open(unit=63,file='junk.out',status='unknown')

! Initialize total helium array to be zero
he_total = 0.0d0


! Calculate helium concentration at each time step by solving implicit finite difference
! approximation to the diffusion equation. The degree of implicitness is beta.
do itime = istart,nresamp

    ! if (printProgress) then
    !     if (mod(itime,nresamp/100).eq.1) then
    !         write(*,1000,advance='no') 'calc_apatite_he_age:',itime*100/nresamp,char(13)
    !     endif
    !     if (itime.eq.nresamp) then
    !         write(*,1000) 'calc_apatite_he_age:',itime*100/nresamp,''
    !     endif
    !     1000 format (1X,A,' progress: [',I3,'% complete]',A)
    ! endif


    !********************!
    ! Produce new helium !
    !********************!
    ! Update the helium production rate (mol/m^3/s)
    call calc_he_production_rate(mol_th232,mol_u235,mol_u238,volume_sphere,production_rate(1))
    production_rate = production_rate(1)
    ! production_rate = 0.0d0

    ! Update the helium concentration at all nodes with newly produced radiogenic helium (mol/m^3)
    ! write(*,*) 'calc_apatite_he_age: producing radiogenic helium'
    do i = 1,nnodes_sphere
        he_conc(i) = he_conc(i) + production_rate(i)*dt_seconds_resamp
    enddo


    !***********************************************!
    ! Check whether diffusion calculation is needed !
    !***********************************************!
    tau = diffusivity(itime)*1.0d0*ma2s/radius_meters**2
    if (tau.gt.tau_max) then
        runDiffusion = .false.
    else
        runDiffusion = .true.
    endif


    !*********!
    ! DIFFUSE !
    !*********!
    if (runDiffusion) then

        !******************************************!
        ! Substitute helium concentration function !
        !******************************************!
        ! Make the substitution for solving differential equation in spherical coordinates
        ! write(*,*) 'calc_apatite_he_age: multiplying concentration by node radii'
        do i = 1,nnodes_sphere
            he_conc_sub(i) = he_conc(i)*radius_shell(i)
            ! print *,i,radius_shell(i),he_conc(i),he_conc_sub(i)
        enddo

        ! Make sure boundary conditions are correct
        he_conc_sub(1) = -he_conc_sub(2)
        he_conc_sub(nnodes_sphere) = he_conc_surf*radius_shell(nnodes_sphere)


        !**************************************************!
        ! Solve for helium concentration at next time step !
        !**************************************************!
        ! Calculate diffusion number at current time step
        diffusion_number = diffusivity(itime)*dt_seconds_resamp/dr_meters**2

        ! write(*,*) 'diffusion_number=',diffusion_number
        call diffusion_step_dirichlet(he_conc_sub,nnodes_sphere,diffusion_number,beta)


        !***************************!
        ! Calculate moles of helium !
        !***************************!
        ! Calculate the new helium concentration at each node ("unsubstitute")
        ! write(*,*) 'calc_apatite_he_age: calculating helium concentration at each node'
        do i = 1,nnodes_sphere
            he_conc(i) = he_conc_sub(i)/radius_shell(i)
            if (he_conc(i).gt.1d20) then
                write(0,*) 'calc_apatite_he_age: looks like you ran into a stability problem...'
                write(0,'(5X,A,I14)')     'number of nodes:        ',nnodes_sphere
                write(0,'(5X,A,1PE14.6)') 'spatial step size (m):  ',dr_meters
                write(0,'(5X,A,1PE14.6)') 'resampled timestep (s): ',dt_seconds_resamp
                write(0,'(5X,A,F14.6)')   'timestep scale factor:  ',dt_var
                write(0,'(5X,A,F14.3)')   'tau | 1 Ma:             ',tau
                write(0,'(5X,A,F14.3)')   'tau_max:                ',tau_max
                write(0,'(5X,A,F14.3)')   'beta:                   ',beta
                write(0,*) 'Options:'
                write(0,*) '    1. Reduce timestep scale factor (-ahe:resamp)'
                write(0,*) '    2. Reduce number of nodes (-ahe:nnodes)'
                write(0,*) '    3. Increase implicitness factor (-ahe:beta)'
                call error_exit(1)
            endif
        enddo

        ! Calculate moles of helium in each shell
        ! write(*,*) 'calc_apatite_he_age: calculating moles of helium in each shell'
        do i = 2,nnodes_sphere-1
            he_mol(i) = volume_shell(i)*he_conc(i)
        enddo

        ! Calculate total moles of helium in the grain
        ! write(*,*) 'calc_apatite_he_age: calculating total helium in the grain'
        do i = 2,nnodes_sphere-1
            he_total(itime) = he_total(itime) + he_mol(i) ! moles
        enddo

    else

        !**************************************!
        ! Set helium to 0 at high temperatures !
        !**************************************!
        he_conc = 0.0d0
        he_total(itime) = 0.0d0

    endif

    !******************************************!
    ! Reduce amount of parent that has decayed !
    !******************************************!
    ! Update amount of radioactive parent elements
    mol_th232 = mol_th232*exp(-decay_th232*dt_seconds_resamp)
    mol_u235 = mol_u235*exp(-decay_u235*dt_seconds_resamp)
    mol_u238 = mol_u238*exp(-decay_u238*dt_seconds_resamp)

enddo


!----
! Calculate (U-Th)/He age of grain
!----
! write(*,*) 'calc_apatite_he_age: calculating (U-Th)/He age'
! print *,'mol_th232',mol_th232
! print *,'mol_u235',mol_u235
! print *,'mol_u238',mol_u238
! print *,'mol_he4',he_total(nresamp)
call calc_u_th_he_age(mol_th232,mol_u235,mol_u238,he_total(nresamp),age)



!----
! Clean up
!----
if (allocated(he_conc)) then
    deallocate(he_conc)
endif
if (allocated(he_total)) then
    deallocate(he_total)
endif



return
end subroutine




end module






