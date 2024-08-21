!--------------------------------------------------------------------------------------------------!
! Module Radiogenic Helium                                                                         !
!                                                                                                  !
! Authors:                                                                                         !
!     - Matt Herman (current Modern Fortran version, i.e., what you are looking at right now!)     !
!     - Rachel Piotraschke (original version; PSU MS thesis)                                       !
!     - Kevin Furlong (original version; supervisor)                                               !
!                                                                                                  !
! This module contains variables and subroutines for the production and diffusion of radiogenic    !
! helium produced from the decay of uranium and thorium, and the analysis of helium-bearing        !
! mineral ages.                                                                                    !
!                                                                                                  !
! References:                                                                                      !
! Farley, K.A. (2000). Helium diffusion from apatite: General behavior as illustrated by Durango   !
!     fluorapatite. Journal of Geophysical Research, 105, B2, 2903-2914.                           !
! Piotraschke, R. (2012). Thermal and Geologic Constraints on the Cretaceous-to-Neogene Tectonic   !
!     Development of the Klamath Mountains, Northern California. M.Sc. Thesis, Pennsylvania State  !
!     University.                                                                                  !
! Prohaska, T., Irrgeher, J., Benefield, J., Böhlke, J.K., Chesson, L. A., Coplen, T.B., Ding, T., !
!     Dunn, P.J.H., Gröning, M., Holden, N.E., Meijer, H.A.J., Moossen, H., Possolo, A.,           !
!     Takahashi, Y., Vogl, J., Walczyk, T., Wang, J., Wieser, M.E., Yoneda, S., Zhu, X.-K., Meija, !
!     J. (2022). Standard atomic weights of the elements 2021 (IUPAC Technical Report). Pure and   !
!     Applied Chemistry, 94(5), 573–600.                                                           !
! Steiger, R.H., Jager, E. (1977). Subcomission on geochronology: Convention on the use of decay   !
!     constants in geo- and cosmochronology. Earth and Planetary Science Letters, 36, 359-362.     !
! Wang, M., Huang, W.J., Kondev, F.G., Audi, G., Naimi, S. (2021). The AME2020 atomic mass         !
!     evaluation. Chinese Physics C, 45(3), 030003.                                                !
!--------------------------------------------------------------------------------------------------!


module radiogenic_helium

    ! Universal gas constant
    double precision, parameter :: univ_gas_constant = 1.9858775d-3             ! [kcal/mol/K]


    ! Helium diffusion in apatite (Farley, 2000: p. 2910)
    ! Note: helium diffusivity in apatite ranges from 0.0008 for c-perpendicular diffusion to
    ! 0.0130 for c-parallel diffusion
    double precision, parameter :: he_diffusivity_apatite = 0.005d0             ! [m^2/s]
    double precision, parameter :: he_activation_energy_apatite = 32.9d0        ! [kcal/mol]


    ! Atomic masses (AME 2020; Wang et al., 2021)
    double precision, parameter :: atomic_mass_th232_ame2020 = 232.0380536d0
    double precision, parameter :: atomic_mass_u235_ame2020  = 235.0439281d0
    double precision, parameter :: atomic_mass_u238_ame2020  = 238.0507869d0

    ! Atomic masses (IUPAC 2021; Prohaska et al., 2022)
    double precision, parameter :: atomic_mass_th232_iupac2021 = 232.03805d0
    double precision, parameter :: atomic_mass_u235_iupac2021  = 235.04393d0
    double precision, parameter :: atomic_mass_u238_iupac2021  = 238.05079d0


    ! Millions of years to seconds
    double precision, parameter :: ma2s = 1.0d6 * 365.25d0 * 24.0d0 * 60.0d0 * 60.0d0

    ! Decay constants (Steiger and Jager, 1977)
    double precision, parameter :: decay_th232 = 4.94750d-5/ma2s                ! [1/s]
    double precision, parameter :: decay_u235  = 9.84850d-4/ma2s                ! [1/s]
    double precision, parameter :: decay_u238  = 1.55125d-4/ma2s                ! [1/s]

    ! Isotopic ratios
    double precision, parameter :: ratio_u238_u235 = 137.88d0


    ! Subroutines
    PUBLIC :: calc_apatite_he_age
    PRIVATE :: calc_he_production_rate
    PRIVATE :: calc_u_th_he_age


contains


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!------------------------------------- He-IN MINERAL ROUTINES -------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!


    subroutine calc_apatite_he_age( &
        n,                          &
        temp_celsius,               &
        dep_km,                     &
        dt_ma,                      &
        dt_var,                     &
        radius_microns,             &
        nnodes,                     &
        beta,                       &
        tau_max,                    &
        verbosity,                  &
        age                         &
    )

    !----
    ! Calculate the concentration of helium in an apatite grain given a temperature-time history.
    ! The grain is assumed to be spherical, and the concentrations of uranium and thorium are
    ! assumed to be homogeneous. The finite difference approximation for the production-diffusion
    ! equation in spherical coordinates is used to determine the helium in the grain.
    !
    ! Inputs:
    !   n:                number of time steps
    !   temp_celsius:     n-point double precision array of temperatures (Celsius)
    !   dep_km:           n-point double precision array of depths, positive up (km)
    !   dt_ma:            time step size for the temperature history (Ma)
    !   dt_var:           variable defining new timestep for diffusion calculation
    !   radius_microns:   radius of the (spherical) apatite grain (microns)
    !   nnodes:           number of spatial nodes + 2 boundary condition nodes
    !   beta:             implicitness weight in finite difference solution (0-1)
    !   tau_max:          maximum dimensionless time to retain He and calculate diffusion (>1 Ma)
    !   verbosity:        level of output to print
    !
    ! Output:
    !   age:              age of grain (Ma)
    !----

    use diffusion, only:              &
        nnodes_sphere,                &
        radius_shell,                 &
        volume_shell,                 &
        volume_sphere,                &
        init_spherical_node_geometry, &
        diffusion_step_dirichlet


    implicit none

    ! Arguments
    integer, intent(in) :: n                                    ! number of input timesteps
    double precision, intent(in) :: temp_celsius(n)             ! temperature history [C]
    double precision, intent(in) :: dep_km(n)                   ! depth history [km]
    double precision, intent(in) :: dt_ma                       ! input timestep size [Ma]
    double precision, intent(in) :: dt_var                      ! >0: dt reduction; <0: timestep [Ma]
    double precision, intent(in) :: radius_microns              ! grain radius [microns]
    integer, intent(in) :: nnodes                               ! number of spatial nodes
    double precision, intent(in) :: beta                        ! implicitness parameter
    double precision, intent(in) :: tau_max                     ! max dimensionless time to keep He
    integer, intent(in) :: verbosity                            ! verbosity level [0-2]
    double precision, intent(out) :: age                        ! (U-Th)/He age [Ma]

    ! Local variables
    integer :: i
    integer :: j
    integer :: itime
    integer :: istart
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
    double precision :: he_conc_surf                            ! He surface BC [mol/m^3]
    double precision :: he_conc_init                            ! He initial BC [mol/m^3]
    double precision, allocatable :: he_conc(:)                 ! He molar concentration [mol/m^3]
    double precision, allocatable :: he_conc_sub(:)             ! radius x He concentration
    double precision, allocatable :: he_mol(:)                  ! moles of He at each node [mol]
    double precision, allocatable :: he_total(:)                ! total moles He [mol]
    double precision :: diffusion_number                        ! diffusivity x dt / dr^2
    double precision :: arg                                     ! dummy variable for exp arguments
    double precision :: tau                                     ! Non-dimensional time
    logical :: runDiffusion                                     ! Diffusion switch



    !----
    ! Set up finite difference spatial grid (assuming spherical apatite grain)
    !----
    if (verbosity.ge.2) then
        write(*,*) 'calc_apatite_he_age: setting up finite difference spatial grid'
    endif

    ! Generate spherical nodes and geometric variables (diffusion module):
    !   - nnodes_sphere           spatial nodes + 2 BC nodes (at sphere center and beyond sphere edge)
    !   - radius_shell [m]        distance of each node from center of sphere [m]
    !   - volume_shell [m^3]      volume of spherical shell for each node [m^3]
    !   - volume_sphere [m^3]     volume of entire sphere [m^3]
    radius_meters = radius_microns/1.0d6              ! Radius of grain [meters]
    dr_microns = radius_microns/dble(nnodes-2)        ! Spatial step size [microns]
    dr_meters = dr_microns*1.0d-6                     ! Spatial step size [meters]
    call init_spherical_node_geometry( &
        radius_meters,                 &
        dr_meters,                     &
        nnodes                         &
    )
    if (verbosity.ge.3) then
        write(*,*) '    radius_meters=',radius_meters
        write(*,*) '    dr_meters=    ',dr_meters
        write(*,*) '    nnodes_sphere=',nnodes_sphere
        write(*,*) '    volume_sphere=',volume_sphere
        write(*,*) '    Node      Shell_Radius      Shell_Volume'
        do i=1,nnodes_sphere
            write(*,'(5X,I4,2(4X,E14.6))') i,radius_shell(i),volume_shell(i)
        enddo
    endif


    !----
    ! Set up finite difference time stepping
    !----
    if (verbosity.ge.2) then
        write(*,*) 'calc_apatite_he_age: resampling temperature history onto new time grid'
    endif

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

    ! Calculate number of resampled time steps
    nresamp = int(dble(n)*dt_seconds/dt_seconds_resamp)

    if (verbosity.ge.3) then
        write(*,*) '    dt_ma=       ',dt_ma
        write(*,*) '    dt_seconds=  ',dt_seconds
        write(*,*) '    dt_ma_resamp=',dt_seconds_resamp
        write(*,*) '    n=           ',n
        write(*,*) '    nresamp=     ',nresamp
    endif


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
    if (verbosity.ge.2) then
        write(*,*) 'calc_apatite_he_age: calculating helium diffusivity over temperature history'
    endif

    ! Diffusivity at each time step (Farley, 2000) [m^2/s]
    if (.not.allocated(diffusivity)) then
        allocate(diffusivity(nresamp))
    endif
    diffusivity = 0.0d0
    do i = 1,nresamp
        arg = -he_activation_energy_apatite/(univ_gas_constant*temp_kelvin(i))
        diffusivity(i) = he_diffusivity_apatite * exp(arg)
    enddo
    if (verbosity.ge.3) then
        write(*,*) '    min(diffusivity)=     ',minval(diffusivity)
        write(*,*) '    max(diffusivity)=     ',maxval(diffusivity)
        write(*,*) '    min(diffusion_number)=',minval(diffusivity)*dt_seconds_resamp/dr_meters**2
        write(*,*) '    max(diffusion_number)=',maxval(diffusivity)*dt_seconds_resamp/dr_meters**2
    endif

    !----
    ! Initialize radiogenic parents and helium production rate array
    !----
    if (verbosity.ge.2) then
        write(*,*) 'calc_apatite_he_age: initializing radiogenic parents and He production rate array'
    endif

    ! Set moles of uranium and thorium, assuming uniform distribution
    mass_th232 = 0.1d-9                               ! Initial mass of thorium-232 (grams)
    mass_u238 = 0.1d-9                                ! Initial mass of uranium-238 (grams)
    mol_th232 = mass_th232/atomic_mass_th232_ame2020  ! Initial moles of thorium-232
    mol_u238 = mass_u238/atomic_mass_u238_ame2020     ! Initial moles of uranium-238
    mol_u235 = mol_u238/ratio_u238_u235               ! Initial moles of uranium-235

    ! Initialize helium production rate array (production rate recalculated in each timestep)
    if(.not.allocated(production_rate)) then
        allocate(production_rate(nnodes_sphere))
    endif
    production_rate = 0.0d0



    !----
    ! Initialize helium concentration arrays
    !----
    if (verbosity.ge.2) then
        write(*,*) 'calc_apatite_he_age: initializing helium concentration'
    endif

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
    if (verbosity.ge.3) then
        write(*,*) '    he_conc_init=',he_conc_init
        write(*,*) '    he_conc_surf=',he_conc_surf
    endif




    !----
    ! Determine helium concentration in the apatite grain over temperature history
    !----
    if (verbosity.ge.2) then
        write(*,*) 'calc_apatite_he_age: determining helium concentration over time'
    endif
    ! open(unit=63,file='junk.out',status='unknown')

    ! Initialize total helium over time [mol] to be zero
    he_total = 0.0d0


    ! Calculate helium concentration at each time step by solving implicit finite difference
    ! approximation to the diffusion equation. The degree of implicitness is beta.
    do itime = istart,nresamp

        if (verbosity.ge.1) then
            if (mod(itime,nresamp/100).eq.1) then
                write(*,1000,advance='no') 'calc_apatite_he_age:',itime*100/nresamp,char(13)
            endif
            if (itime.eq.nresamp) then
                write(*,1000) 'calc_apatite_he_age:',itime*100/nresamp,''
            endif
            1000 format (1X,A,' progress: [',I3,'% complete]',A)
        endif


        !********************!
        ! Produce new helium !
        !********************!
        ! Update the helium production rate [mol/m^3/s]
        call calc_he_production_rate( &
            mol_th232,                &
            mol_u235,                 &
            mol_u238,                 &
            volume_sphere,            &
            production_rate(1)        &
        )
        production_rate = production_rate(1)
        ! production_rate = 0.0d0

        ! Update the helium concentration at all nodes with new radiogenic helium [mol/m^3]
        ! write(*,*) 'calc_apatite_he_age: producing radiogenic helium'
        do i = 1,nnodes_sphere
            he_conc(i) = he_conc(i) + production_rate(i)*dt_seconds_resamp
        enddo


        !***********************************************!
        ! Check whether diffusion calculation is needed !
        !***********************************************!
        tau = diffusivity(itime)*ma2s/radius_meters**2
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
            call diffusion_step_dirichlet( &
                he_conc_sub,               &
                nnodes_sphere,             &
                diffusion_number,          &
                beta                       &
            )


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
                he_mol(i) = volume_shell(i)*he_conc(i)        ! moles
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
        mol_u235  = mol_u235* exp(-decay_u235* dt_seconds_resamp)
        mol_u238  = mol_u238* exp(-decay_u238* dt_seconds_resamp)

    enddo



    !----
    ! Calculate (U-Th)/He age of grain
    !----
    if (verbosity.ge.2) then
        write(*,*) 'calc_apatite_he_age: calculating (U-Th)/He age'
    endif
    if (verbosity.ge.3) then
        write(*,*) '    mol_th232',mol_th232
        write(*,*) '    mol_u235 ',mol_u235
        write(*,*) '    mol_u238 ',mol_u238
        write(*,*) '    mol_he4  ',he_total(nresamp)
    endif

    ! Solve radiometric age equation
    call calc_u_th_he_age( &
        mol_th232,         &
        mol_u235,          &
        mol_u238,          &
        he_total(nresamp), &
        age                &
    )
    if (verbosity.ge.1) then
        write(*,*) '    age',age
    endif



    !----
    ! Clean up
    !----
    if (allocated(temp_kelvin)) then
        deallocate(temp_kelvin)
    endif
    if (allocated(diffusivity)) then
        deallocate(diffusivity)
    endif
    if (allocated(production_rate)) then
        deallocate(production_rate)
    endif
    if (allocated(he_conc)) then
        deallocate(he_conc)
    endif
    if (allocated(he_conc_sub)) then
        deallocate(he_conc_sub)
    endif
    if (allocated(he_mol)) then
        deallocate(he_mol)
    endif
    if (allocated(he_total)) then
        deallocate(he_total)
    endif


    return
    end subroutine




!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!----------------------------- (U-Th)/He SYSTEM SUBROUTINES ---------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!



    subroutine calc_he_production_rate( &
        mol_th232,                      &
        mol_u235,                       &
        mol_u238,                       &
        volume,                         &
        p                               &
    )
    !----
    ! Determine the radiogenic helium production rate [mol/m^3/s] in a volume, given the amount of
    ! radioactive thorium-232, uranium-235, and uranium-238 [moles]. The distribution of parent
    ! nuclides within the volume is assumed to be homogeneous.
    !
    ! Inputs:
    !   mol_th232:  moles of thorium-232
    !   mol_u235:   moles of uranium-235
    !   mol_u238:   moles of uranium-238
    !   volume:     volume [m^3]
    !
    ! Outputs:
    !   p:          helium production rate [mol/m^3/s]
    !----

    implicit none

    ! Arguments
    double precision, intent(in) :: mol_th232
    double precision, intent(in) :: mol_u235
    double precision, intent(in) :: mol_u238
    double precision, intent(in) :: volume
    double precision, intent(out) :: p

    ! Local variables
    double precision :: mol_conc_th232
    double precision :: mol_conc_u235
    double precision :: mol_conc_u238


    ! Molar concentration (mol/m^3)
    mol_conc_th232 = mol_th232/volume
    mol_conc_u235 = mol_u235/volume
    mol_conc_u238 = mol_u238/volume

    ! Helium production rate (mol/m^3/s)
    ! thorium-232 -> lead-208: 6 alpha particles (helium nuclei)
    ! uranium-235 -> lead-207: 7 alpha particles
    ! uranium-238 -> lead-206: 8 alpha particles
    p = 8.0d0 * mol_conc_u238 * decay_u238 &
            + 7.0d0 * mol_conc_u235 * decay_u235 &
            + 6.0d0 * mol_conc_th232 * decay_th232


    return
    end subroutine



!--------------------------------------------------------------------------------------------------!



    subroutine calc_u_th_he_age( &
        mol_th232,               &
        mol_u235,                &
        mol_u238,                &
        mol_he4,                 &
        age                      &
    )
    !----
    ! Determine the (U-Th)/He age (in Ma) based on the amount of parent uranium and thorium and
    ! daughter helium present in a sample by solving the radiometric age equation for the (U-Th)/He
    ! system. This equation does not have a closed form solution, so here it is solved by
    ! searching for a zero crossing with progressively shorter intervals.
    !
    ! Inputs:
    !   mol_th232:  moles of thorium-232
    !   mol_u235:   moles of uranium-235
    !   mol_u238:   moles of uranium-238
    !   mol_he4:    moles of helium-4
    !
    ! Outputs:
    !   age:        Age [Ma]
    !----

    implicit none

    ! Arguments
    double precision, intent(in) :: mol_th232
    double precision, intent(in) :: mol_u235
    double precision, intent(in) :: mol_u238
    double precision, intent(in) :: mol_he4
    double precision, intent(out) :: age

    ! Local variables
    double precision :: t_ma
    double precision :: dt_ma
    double precision :: tmax_ma
    double precision :: t_seconds
    double precision :: dt_seconds
    double precision :: tmax_seconds
    double precision :: f


    ! Initialize age zero-crossing search parameters
    t_ma = 0.0d0              ! Initial age [Ma]
    tmax_ma = 4.5d3           ! Maximum age [Ma]
    dt_ma = 1.0d0             ! Starting age step size [Ma]

    ! Convert to seconds
    t_seconds = t_ma*ma2s
    tmax_seconds = tmax_ma*ma2s
    dt_seconds = dt_ma*ma2s

    ! Solve for age by searching for zeros of helium production function
    do while (t_seconds.le.tmax_seconds)

        ! Age function, f, has the following characteristics:
        !     (a) starts negative
        !     (b) is monotonically increasing
        !     (c) its zero is the age
        f = 8.0d0*mol_u238*(exp(decay_u238*t_seconds)-1.0d0) + &
            7.0d0*mol_u235*(exp(decay_u235*t_seconds)-1.0d0) + &
            6.0d0*mol_th232*(exp(decay_th232*t_seconds)-1.0d0) - &
            mol_he4

        ! Did the age function, f, cross zero?
        if (f.gt.0.0d0) then
            t_seconds = t_seconds - dt_seconds ! Step back to the last age before crossing zero
            dt_seconds = dt_seconds/10.0d0     ! Decrease the time step size
            if (dt_seconds.lt.1d-6*ma2s) then  ! Finish when time step is less than 1 year
                exit
            endif
        endif

        ! Advance to next time
        t_seconds = t_seconds + dt_seconds

    enddo

    ! Calculate age in Ma
    age = t_seconds/ma2s

    return
    end subroutine


end module