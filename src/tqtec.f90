!==================================================================================================!
!                                                                                                  !
! TQ[Tec|Clim]                                                                                     !
!                                                                                                  !
! TQTec (Temperature, Heat Flow, Tectonics) and TQClim (Temperature, Heat Flow, Climate)           !
! calculate the one-dimensional transient thermal field within an area that experiences geologic,  !
! tectonic, and/or climatic events. Processes that can be modeled include:                         !
!     - Burial                                                                                     !
!     - Uplift/Erosion                                                                             !
!     - Thrust Faulting (monitoring either the hanging wall or footwall)                           !
!     - Bulk Thickening                                                                            !
!     - Surface Temperature Variations                                                             !
!                                                                                                  !
! Authors:                                                                                         !
!     - Kevin Furlong (original Fortran 77 program)                                                !
!     - Matt Herman (Modern Fortran version, i.e., what you are looking at right now!)             !
!     - Chris Guzofski, Matt Legg, Rachel Piotraschke (PSU MS theses)                              !
!                                                                                                  !
!==================================================================================================!



module tqtec


! Inputs/outputs
character(len=512) :: input_file                  ! name of input file
character(len=8) :: input_mode                    ! how to read input parameters (user, file)
character(len=512) :: output_file                 ! name of output file
character(len=512) :: geotherm_file               ! name of output geotherm file
character(len=512) :: timing_file                 ! name of output tectonic action timing file
integer :: verbosity                              ! level of information to print during execution


! Finite difference parameters
integer :: nnodes                                 ! number of spatial nodes
integer :: nt_total                               ! number of time steps
integer :: istep                                  ! current time step
double precision :: dz                            ! node spacing (tqtec:km, tqclim:m)
double precision :: dt                            ! time step interval (tqtec:Ma, tqclim:yr)
double precision :: r1                            ! finite difference time factor


! Timing
double precision :: t_total                       ! total model time (tqtec:Ma, tqclim:yr)
double precision :: t_geotherm_output             ! time per geotherm output (tqtec:Ma, tqclim:yr)
integer :: nt_geotherm_output                     ! time steps between geotherm outputs


! Nodal parameters
double precision, allocatable :: conductivity(:)  ! conductivity (W/m/K)
double precision, allocatable :: temp(:)          ! temperature (C)
double precision, allocatable :: hp(:)            ! heat production (uW/m^3)


! Time-varying parameters
double precision, allocatable :: hf(:)            ! surface heat flow (mW/m^2)


! Horizons of interest
integer :: nhorizons                              ! number of horizons
double precision, allocatable :: depth(:)         ! depth of horizons (tqtec:km, tqclim:m)
integer, allocatable :: depth_node(:)             ! horizon nodes


! Material properties
integer :: nlayers                                ! number of distinct material layers
double precision, allocatable :: layer(:,:)       ! layer(:,1): depth to top (km)
                                                  ! layer(:,2): thickness (km)
                                                  ! layer(:,3): conductivity (W/m/K)
double precision :: diffusivity                   ! diffusivity (currently a hard-coded constant)
double precision :: cond_base                     ! basal conductivity (W/m/K)


! Boundary conditions
double precision :: temp_surf                     ! surface temperature (C)
double precision :: hp_surf                       ! surface heat production (uW/m^3)
double precision :: hp_dep                        ! depth of heat production (km)
double precision :: hf_surf                       ! surface heat flow (mW/m^2)
double precision :: hf_base                       ! basal heat flow (mW/m^2)
double precision :: dtemp_wo_hp                   ! temp change w/o heat prod (C)
double precision :: temp_factor                   ! temp scaling factor
double precision :: temp_base_adj                 ! temp at node nnodes+1 (c)


! Thermal/tectonic/climate events
logical :: isActionDefined
integer, allocatable :: action_array(:)           ! Tectonic/thermal/climatic actions encoded as
                                                  ! base 2 digits (saved as base 10):
                                                  ! 0|1 0|1 0|1 0|1 0|1
                                                  !  ^   ^   ^   ^   ^--- bury
                                                  !  |   |   |   |------- uplift/erode
                                                  !  |   |   |----------- thrust
                                                  !  |   |--------------- thicken
                                                  !  |------------------- thin

integer :: nburial !------------------------------! number of burial events
double precision, allocatable :: burial_dat(:,:)  ! burial_dat(:,1): start (Ma)
                                                  ! burial_dat(:,2): duration (Ma)
                                                  ! burial_dat(:,3): thickness (km)
                                                  ! burial_dat(:,4): conductivity (W/(m*K))
double precision, allocatable :: burial_cond(:)   ! conductivity of new material at each timestep


integer :: nuplift !------------------------------! number of uplift/erosion events
double precision, allocatable :: uplift_dat(:,:)  ! uplift_dat(:,1): start (Ma)
                                                  ! uplift_dat(:,2): duration (Ma)
                                                  ! uplift_dat(:,3): thickness (km)


integer :: nthrust !------------------------------! number of thrusting events
integer, allocatable :: thrust_step(:)            ! thrusting timestep array
double precision, allocatable :: thrust_dat(:,:)  ! thrust_dat(:,1): start (Ma)
                                                  ! thrust_dat(:,2): upper (1) or lower (2) plate
                                                  ! thrust_dat(:,3): initial base (km)
                                                  ! thrust_dat(:,4): initial emplacement depth (km)
                                                  ! thrust_dat(:,5): initial thickness (km)


integer :: nhfvars !------------------------------! number of surface heat flow variations
double precision, allocatable :: hfvar(:,:)       ! (1) start (2) new heat flow
double precision, allocatable :: bas_grad(:)      ! temperature gradient at the base of the model


integer :: nthicken !-----------------------------! number of thickening/thinning events
double precision, allocatable :: thicken_dat(:,:) ! thicken_dat(:,1): start (Ma)
                                                  ! thicken_dat(:,2): duration (Ma)
                                                  ! thicken_dat(:,3): total thickening amount (km)
                                                  ! thicken_dat(:,4): top of crust (km)
                                                  ! thicken_dat(:,5): initial crustal thickness (km)
integer, allocatable :: thicken_start(:)          ! thicken_start(:): timestep thickening starts
integer, allocatable :: crust_dat(:,:)            ! crust_dat(:,1): top of crust (node)
                                                  ! crust_dat(:,2): bottom of crust (node)
integer:: crust_top
integer:: crust_bot
logical :: thickenHorizons                        ! T: move tracked horizon depths with thickening
integer, allocatable :: horizon_shift(:,:)


integer :: ntempsteps !---------------------------! number of surface temperature step variations
integer :: ntempramps !---------------------------! number of surface temperature ramp variations
integer :: ntempcycles !--------------------------! number of surface temperature ramp variations
double precision, allocatable :: temp_step_dat(:,:)   ! (1) start (2) new temp
double precision, allocatable :: temp_ramp_dat(:,:)   ! (1) start (2) duration (3) temp change
double precision, allocatable :: temp_cycle_dat(:,:)  ! (1) start (2) duration (3) amp (4) freq
double precision, allocatable :: temp_surf_var(:)     ! surface temperature variations over time


! Results array
double precision, allocatable :: results(:,:,:)   ! temperature and depth for each timstep/depth


end module tqtec







!==================================================================================================!
!==================================================================================================!
!==================================================================================================!
!==================================================================================================!
!==================================================================================================!




program main

!----
! Solve for the 1-D transient thermal field defined by boundary conditions, material properties,
! and tectonic/geologic/climate actions
!----


use tqtec


implicit none


! Local variables
integer :: i
integer :: j
integer :: np
integer :: ierr
integer :: iplate
double precision :: cond_surf
character(len=8) :: exec_name
character(len=2) :: time_unit



! Set name of executable and associated variables
#ifdef COMPILE_TQTEC
    exec_name = 'tqtec'
    time_unit = 'Ma'
#elif COMPILE_TQCLIM
    exec_name = 'tqclim'
    time_unit = 'yr'
#else
    exec_name = 'NAME_ERR'
#endif




! Initialize default model parameters
call initialize_defaults()                            ! tqtec_io.f90



! Parse command line options
call gcmdln()                                         ! tqtec_io.f90



if (verbosity.ge.1) then
    write(*,*) trim(exec_name)//': starting'
endif



! Read control file or user input from standard input
call read_model_parameters()                          ! tqtec_io.f90





! Check that model is deep enough to include all depth horizons of interest
if (dz*dble(nnodes).lt.maxval(depth)) then
    write(0,*) trim(exec_name)//': model extent is shallower than deepest horizon'
    call error_exit(1)
endif



! Initialize the temperature, heat flow, and heat production at each node
call initialize_thermal_parameters()                  ! tqtec.f90



! Set up action timing arrays
call setup_action_arrays()                            ! tqtec_actions.f90
call check_actions()                                  ! tqtec_actions.f90



! Print model parameters to standard output
if (verbosity.ge.1) then
    call print_model_parameters()                     ! tqtec_io.f90
endif



! Print the geotherm to a file (-geotherm flag)
if (geotherm_file.ne.'') then                         ! tqtec_io.f90
    call print_geotherm( &
        geotherm_file,   &
        0,               &
        0.0d0,           &
        t_total,         &
        nnodes,          &
        temp,            &
        dz               &
    )
endif



! Step through time and run the finite difference procedure to calculate the temperature at each
! node for each time step

! Initialize timestep
istep = 0

! Until model reaches the final step...
do while (istep.lt.nt_total)


    ! Update the adjusted temperature at the base of the model
    temp_base_adj = temp(nnodes) + bas_grad(istep+1)


    ! Calculate the updated temperatures at each node (***THE MAIN FINITE DIFFERENCE PROCEDURE***)
    call update_temps(nnodes,ierr)                    ! tqtec.f90
    if (ierr.ne.0) then
        write(0,*) trim(exec_name)//': error in update_temps() TRID algorithm at step',istep
        call error_exit(1)
    endif


    ! Increment the time step
    istep = istep + 1
    if (verbosity.le.2) then
        if (istep.lt.nt_total) then
            write(*,'(X,A,I6,A,I6,A)',advance='no') trim(exec_name)//': working on step',istep,&
                                                    ' of',nt_total,char(13)
        else
            write(*,'(X,A,I6,A,I6)') trim(exec_name)//': working on step',istep,' of',nt_total
        endif
    elseif (verbosity.ge.3) then
        if (istep.lt.nt_total.and.action_array(istep).eq.0) then
            write(*,'(X,A,I6,A,I6,A)',advance='no') trim(exec_name)//': working on step',istep,&
                                                    ' of',nt_total,char(13)
        else
            write(*,'(X,A,I6,A,I6)') trim(exec_name)//': working on step',istep,' of',nt_total
        endif
    endif



    ! Thermal/tectonic/climate actions                ! tqtec_actions.f90

    ! Add material to the top of the model (burial)
    if (mod(action_array(istep)/1,2).eq.1) then
        call bury()
    endif

    ! Remove material from the top of the model (uplift + erosion)
    if (mod(action_array(istep)/2,2).eq.1) then
        call erode()
    endif

    ! Duplicate and emplace material (thrust sheet)
    if (mod(action_array(istep)/4,2).eq.1) then
        ! Get thrust event number
        j = 0
        do i = 1,nthrust
            if (thrust_step(i).eq.istep) then
                j = i
                exit
            endif
        enddo
        if (j.eq.0) then
            write(0,*) trim(exec_name)//': error setting thrust number at timestep',istep
        endif
        iplate = int(thrust_dat(j,2))
        if (iplate.eq.0) then
            write(0,*) trim(exec_name)//': error setting thrust number at timestep',istep
        endif
        ! Duplicate material and insert it (thrusting)
        if (iplate.eq.1) then
            ! Track horizons in the hanging wall (upper plate) of the thrust sheet
            call thrust('upper')
        elseif (iplate.eq.2) then
            ! Track horizons in the footwall (lower plate) of the thrust sheet
            call thrust('lower')
        endif
    endif

    ! Add or remove material within model (crustal thickening/thinning)
    if (mod(action_array(istep)/8,2).eq.1.or.mod(action_array(istep)/16,2).eq.1) then
        ! Check to see if this is a new crustal thickening/thinning event
        do i = 1,nthicken
            if (thicken_start(i).eq.istep) then
                ! Reset initial crustal top and bottom nodes
                crust_top = crust_dat(i,1)
                crust_bot = crust_dat(i,2)
                exit
            endif
        enddo
        ! Add or remove a node within the model interior (thickening/thinning)
        if (mod(action_array(istep)/8,2).eq.1) then
            call thicken()
        elseif (mod(action_array(istep)/16,2).eq.1) then
            call thin()
        else
            write(0,*) trim(exec_name)//': somehow you got in the thickening/thinning action block'
            write(0,*) 'But action code',action_array(istep),' does not indicate either action...'
            call error_exit(1)
        endif

    endif



    ! Check if horizons should be shifted at this timestep
    ! NOTE: I want to move this to thicken()/thin(), but when horizon nodes move during a thickening
    ! or thinning event is different from when a node i added to the crust during that event
    if (nthicken.gt.0.and.thickenHorizons) then
        depth_node = depth_node + horizon_shift(:,istep)
    endif


    ! Calculate surface heat flow for this time step
    hf(istep) = (temp(10)-temp(5))/(5.0d0*dz)   ! Temperature gradient near surface
    cond_surf = 0.0d0                           ! Average surface conductivity
    do i = 1,5
        cond_surf = cond_surf + conductivity(i+4)
    enddo
    cond_surf = cond_surf/5.0d0
    hf(istep) = hf(istep)*cond_surf             ! Heat flow = dT/dz * conductivity


    ! Update surface temperature
    temp_surf = temp_surf_var(istep)


    ! Save depths and temperatures of tracked horizons in results array
    do i = 1,nhorizons
        np = depth_node(i)
        if (np.eq.0) then
            results(istep,1,i) = temp_surf
        elseif (np.lt.0) then
            results(istep,1,i) = 0.0d0
        elseif (np.gt.0) then
            results(istep,1,i) = temp(np)
        endif
        results(istep,2,i) = dble(depth_node(i))
    enddo


    ! Print geotherm every nt_geotherm_output steps
    if (geotherm_file.ne.'') then
        if (mod(istep,nt_geotherm_output).eq.0) then
            call print_geotherm( &
                geotherm_file,   &
                istep,           &
                dble(istep)*dt,  &
                t_total,         &
                nnodes,          &
                temp,            &
                dz               &
            )
        endif
    endif



enddo   ! end of timestep loop



! Close the geotherm file if needed
if (geotherm_file.ne.'') then
    close(12)
endif



! Print the results to the output file
call output()                                         ! tqtec_io.f90



if (verbosity.ge.1) then
    write(*,*) trim(exec_name)//': finished'
    write(*,*) 'Results can be found in ',trim(output_file)
    if (geotherm_file.ne.'') then
        write(*,*) 'Geotherms every ',t_geotherm_output,trim(time_unit), &
                   ' can be found in ',trim(geotherm_file)
    endif
endif






end program






















!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!---------------------------------- MODEL PREPARATION ---------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!




subroutine initialize_thermal_parameters()
!----
! Calculate the steady state temperature at each node based on the surface heat flow, surface
! temperature, conductivity, and heat production
!----

use tqtec

implicit none

! Local variables
integer :: i, j
integer :: ntop, nbot
double precision :: hfhp



if (verbosity.ge.2) then
    write(*,*) 'initialize_thermal_parameters: starting'
endif


! Calculate model parameters
nt_total = int(t_total/dt)                            ! number of timesteps
nt_geotherm_output = int(t_geotherm_output/dt)        ! number of timesteps between geotherm outputs
temp_factor = diffusivity*dt/cond_base                ! temperature scaling factor
r1 = diffusivity*dt/(dz*dz)                           ! finite difference time factor
dtemp_wo_hp = (hf_surf - hp_surf*hp_dep)*dz/cond_base ! temperature difference without heat production



! Allocate arrays
allocate(depth_node(nhorizons))                       ! node where each horizon of interest sits
allocate(conductivity(nnodes))                        ! conductivity at each node
allocate(temp(nnodes))                                ! temperature at each node
allocate(hp(nnodes))                                  ! heat production at each node
allocate(hf(nt_total))                                ! surface heat flow at each timestep
allocate(results(nt_total,2,nhorizons))               ! temperature/depth at each time step & horizon



! Initialize the conductivity at each node to be the basal conductivity
do i = 1,nnodes
    conductivity(i) = cond_base
enddo


! If there are multiple layers with different conductivities, then locate each node within a layer
! and set the nodal conductivity to be the corresponding layer conductivity
if (nlayers.gt.0) then
    do i = 1,nlayers
        ! Find node numbers at the top and bottom of the layer
        ntop = int(layer(i,1)/dz)
        nbot = int((layer(i,1)+layer(i,2))/dz)
        ! Set the conductivity at all of the nodes within the layer
        if (ntop+1.eq.nbot) then
            conductivity(ntop+1) = layer(i,3)
        else
            do j = ntop+1,nbot
                conductivity(j) = layer(i,3)
            enddo
        endif
    enddo
endif


! Divide horizon depths by vertical node spacing to place horizons at a node
do i = 1,nhorizons
    depth_node(i) = int(depth(i)/dz)
enddo


! Calculate volumetric heat production at each node based on exponentially decaying heat production
if (hp_dep.gt.0.0d0) then
    do i = 1,nnodes
        hp(i) = hp_surf*exp(-dble(i)*dz/hp_dep)
    enddo
else
    hp = 0.0d0
endif


! Calculate steady state temperature at each node
! Start with node 1 and work from the surface downward
temp(1) = temp_surf + hf_surf*dz/conductivity(1) - hp(1)*dz**2/(2.0d0*conductivity(1))
! Subtract the heat production between nodes to get the heat flow at node 2
hfhp = hp(1)*dz
hf_base = hf_surf - hfhp
! Work downward for all nodes
do i = 2,nnodes
    temp(i) = temp(i-1) + hf_base*dz/conductivity(i) - hp(i)*dz**2/(2.0d0*conductivity(i))
    hfhp = hp(i)*dz
    hf_base = hf_base - hfhp
enddo


! Initialize surface heat flow at time step 1
hf(1) = (conductivity(1)*(temp(1)-temp_surf))/dz


if (verbosity.ge.2) then
    write(*,*) 'initialize_thermal_parameters: finished'
endif



return


end subroutine














!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!------------------------------- FINITE DIFFERENCE PROCEDURE --------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!



subroutine update_temps(nnodes,ierr)
!----
! A combination of old subroutines MAT and TRID
!
! This procedure solves the finite difference approximation to the 1-D heat equation:
!
!                dT              d             dT
!    rho * Cp * ----  =  q  -  ---- [  k(z) * ---- ]
!                dt             dz             dz
!
!    where:
!        T: temperature
!        t: time
!        z: depth
!        q: heat production
!        rho: density
!        Cp: heat capacity
!        k(z): conductivity
!
!----

use tqtec, only: verbosity, &
                 r1, &
                 conductivity, &
                 temp, &
                 hp, &
                 temp_surf, &
                 temp_factor, &
                 temp_base_adj, &
                 dt, &
                 diffusivity

implicit none

! Arguments
integer :: nnodes
integer :: ierr

! Local variables
integer :: i, k, l
double precision :: bet(nnodes), gam(nnodes)
double precision :: c(nnodes), d(nnodes), e(nnodes)
double precision :: a(nnodes,3)
double precision :: temp_new(nnodes)
double precision :: tmp



if (verbosity.ge.3) then
    write(*,*) '    update_temps: starting'
endif



! Calculate nodal coefficients corresponding to variations in conductivity
gam(1) = 1.0d0/conductivity(1)
bet(1) = 1.0d0/(gam(1)+gam(1))
do i = 2,nnodes
    gam(i) = 1.0d0/conductivity(i)
    bet(i) = 1.0d0/(gam(i)+gam(i-1))
enddo


! Calculate finite difference terms used to determine the temperatures at the current time step
! and in the tridiagonal matrix of the finite difference scheme
c(1) = -r1*bet(1)*gam(1)
e(1) = -r1*bet(2)*gam(1)
d(1) = 1.0d0 - c(1) - e(1)
a(1,1) = -c(1)
a(1,2) = 2.0d0 - d(1)
a(1,3) = -e(1)
do i = 2,nnodes-1
    c(i) = -r1*bet(i)*gam(i)
    e(i) = -r1*bet(i+1)*gam(i)
    d(i) = 1.0d0 - c(i) - e(i)
    a(i,1) = -c(i)
    a(i,2) = 2.0d0 - d(i)
    a(i,3) = -e(i)
enddo
c(nnodes) = -r1*bet(nnodes)*gam(nnodes)
e(nnodes) = c(nnodes)
d(nnodes) = 1.0d0 - 2.0d0*c(nnodes)
a(nnodes,1) = -c(nnodes)
a(nnodes,2) = 2.0d0 - d(nnodes)
a(nnodes,3) = -e(nnodes)


! Calculate temperatures at each node for the current time step in the finite difference scheme
temp_new(1) = a(1,2)*temp(1) + a(1,3)*temp(2) + 2.0d0*a(1,1)*temp_surf + &
              diffusivity*dt*hp(1)/conductivity(1)
do i = 2,nnodes-1
    temp_new(i) = a(i,1)*temp(i-1) + a(i,2)*temp(i) + a(i,3)*temp(i+1) + &
                  diffusivity*dt*hp(i)/conductivity(i)
enddo
temp_new(nnodes) = a(nnodes,1)*temp(nnodes-1) + a(nnodes,2)*temp(nnodes) + &
               2.0d0*a(nnodes,3)*temp_base_adj + temp_factor*hp(nnodes)

! Update temperatures of current time step in global temperature array
temp = temp_new



! END SUBROUTINE MAT
! BEGIN SUBROUTINE TRID



! Initialize error flag
ierr = 0


!----
! The c, d, and e arrays are the components of the tridiagonal matrix for calculating the
! temperatures at the next time step. The matrix (A) is:
!
! [ d(1) e(1)  0    0    ...    0      0     0     ]
! [ c(2) d(2) e(2)  0           0      0     0     ]
! [  0   c(3) d(3) e(3)         0      0     0     ]
! [  :                    .                  :     ]
! [  0    0    0    0         c(n-1) d(n-1) e(n-1) ]
! [  0    0    0    0    ...    0    c(n)   d(n)   ]
!
! Then, A * temp(next_step) = temp(current_step)
!----

! Solve matrix equation for temperature at next time step
! Update c array for node 1
c(1) = d(1)
if (nnodes-1.ge.1) then
    ! Update d and e arrays for node 1 and last node
    d(1) = e(1)
    e(1) = 0.0d0
    e(nnodes) = e(1)

    ! Loop through nodes and do forward substitution for matrix solution
    do k = 1,nnodes-1
        ! Flip equation order?
        if (abs(c(k+1)).ge.abs(c(k))) then
            tmp = c(k+1)
            c(k+1) = c(k)
            c(k) = tmp
            tmp = d(k+1)
            d(k+1) = d(k)
            d(k) = tmp
            tmp = e(k+1)
            e(k+1) = e(k)
            e(k) = tmp
            tmp = temp(k+1)
            temp(k+1) = temp(k)
            temp(k) = tmp
        endif
        ! Problem solving matrix equation if c is 0
        if (abs(c(k)).lt.1.0d-8) then
            ierr = k
            return
        endif
        ! Decomposition and forward substitution
        tmp = -c(k+1)/c(k)
        c(k+1) = d(k+1) + tmp*d(k)
        d(k+1) = e(k+1) + tmp*e(k)
        e(k+1) = 0.0d0
        temp(k+1) = temp(k+1) + tmp*temp(k)
    enddo
endif

! Again, c should not be 0
if (abs(c(nnodes)).lt.1.0d-8) then
    ierr = nnodes
    return
endif

! Update last node
temp(nnodes) = temp(nnodes)/c(nnodes)

! Not much else to do with only one node
if (nnodes.eq.1) then
    return
endif

! Update second to last node
temp(nnodes-1) = (temp(nnodes-1)-d(nnodes-1)*temp(nnodes))/c(nnodes-1)

! Not much else to do with only two nodes
if (nnodes-2.lt.1) then
    return
endif


! Backsubstitution to calculate temperatures at next step
do l = 1,nnodes-2
    k = nnodes-2-l+1
    temp(k) = (temp(k)-d(k)*temp(k+1)-e(k)*temp(k+2))/c(k)
enddo


if (verbosity.ge.3) then
    write(*,*) '    update_temps: finished'
endif


return
end subroutine



