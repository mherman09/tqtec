!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------- PREPARE ACTION ARRAYS ------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!


subroutine setup_action_arrays()
!----
! Construct arrays to control tectonic action events
!----


use tqtec


implicit none


! Local variables
integer :: i, ct
integer :: ihor
integer :: ndep
integer :: j, jbeg, jend
integer :: itop
integer :: ithick
integer :: nstart, nduration, nthick !, top_crust, crust_thick, bot_crust
double precision :: rate, arg
double precision :: hf_surf_var(nt_total)
double precision :: delta_depth
integer :: intqt(nhfvars)
integer, allocatable :: total_shift(:)
double precision :: thickening_ratio
double precision :: amp
double precision :: freq
double precision, parameter :: pi = 4.0d0*atan(1.0d0)
logical :: isThinning
character(len=6) :: exec_name



! Set executable name
#ifdef COMPILE_TQTEC
exec_name = 'tqtec'
#elif COMPILE_TQCLIM
exec_name = 'tqclim'
#endif



if (verbosity.ge.2) then
    write(*,*) 'setup_action_arrays: starting'
endif



! Initialize global action arrays
call init_action_arrays()



! Initialize local action arrays
if (allocated(total_shift)) then
    deallocate(total_shift)
endif
if (nthicken.gt.0.and.thickenHorizons) then
    allocate(total_shift(nhorizons))
endif
if (nthicken.gt.0.and.thickenHorizons) then
    total_shift = 0
endif





!**************************************************************************************************!
! *** Burial Events *******************************************************************************!
!**************************************************************************************************!
do i = 1,nburial
    nstart = int(burial_dat(i,1)/dt)              ! Starting timestep
    nduration = int(burial_dat(i,2)/dt)           ! Duration in timesteps
    nthick = int(burial_dat(i,3)/dz)              ! Thickness in nodes
    rate = dble(nthick)/dble(nduration)           ! Rate in nodes/timestep
    jbeg = nstart + 1                             ! First timestep
    jend = nstart + nduration                     ! Last timestep
    ct = 0                                        ! Initialize counter for burial increments
    do j = jbeg,jend
        arg = dble(j-nstart)*rate - dble(ct)      ! Test for burying at this timestep
        if (arg.ge.1.0d0) then
            action_array(j) = action_array(j) + 1
            burial_cond(j) = burial_dat(i,4)
            ct = ct + 1
        endif
    enddo
enddo



!**************************************************************************************************!
! *** Uplift/Erosion Events ***********************************************************************!
!**************************************************************************************************!
do i = 1,nuplift
    nstart = int(uplift_dat(i,1)/dt)              ! Starting timestep
    nduration = int(uplift_dat(i,2)/dt)           ! Duration in timesteps
    nthick = int(uplift_dat(i,3)/dz)              ! Thickness in nodes
    rate = dble(nthick)/dble(nduration)           ! Rate in nodes/timestep
    jbeg = nstart + 1                             ! First timestep
    jend = nstart + nduration                     ! Last timestep
    if (jbeg.le.0) then
        write(0,*) trim(exec_name)//': uplift event timing error'
        write(0,*) 'Start time must be >= 0'
        call error_exit(1)
    endif
    if (jend.gt.nt_total) then
        write(0,*) trim(exec_name)//': uplift event timing error'
        write(0,*) 'Duration set too long, event runs past end of model'
        call error_exit(1)
    endif
    ct = 0                                        ! Initialize counter for uplift increments
    do j = jbeg,jend
        arg = dble(j-nstart)*rate - dble(ct)      ! Test for eroding at this timestep
        if (arg.ge.1.0d0) then
            action_array(j) = action_array(j) + 2
            ct = ct + 1
        endif
    enddo
enddo



!**************************************************************************************************!
! *** Thrust Events *******************************************************************************!
!**************************************************************************************************!
do i = 1,nthrust
    nstart = int(thrust_dat(i,1)/dt) + 1          ! Timestep of thrust faulting
    action_array(nstart) = action_array(nstart) + 4
    thrust_step(i) = nstart                       ! Save thrust timestep for each thrust action
    ! THTYPE(I): thrust_dat(i,2)
enddo



!**************************************************************************************************!
! *** Surface Heat Flow Changes *******************************************************************!
!**************************************************************************************************!
hf_surf_var = hf_surf                             ! Initialize surface heat flow at initial value
do i = 1,nhfvars
    intqt(i) = int(hfvar(i,1)/dt)                 ! Timestep of heat flow change
enddo
do i = 1,nhfvars-1
    jbeg = intqt(i)
    jend = intqt(i+1)
    do j = jbeg,jend
        hf_surf_var(j) = hfvar(i,2)               ! Update surface heat flow
    enddo
enddo
if (nhfvars.ge.1) then
    jbeg = intqt(nhfvars)
    jend = nt_total
    do j = jbeg,jend
        hf_surf_var(j) = hfvar(nhfvars,2)         ! Update surface heat flow (to end of model)
    enddo
endif

! Check that heat production is never greater than surface heat flow
if (hp_surf*hp_dep.ge.maxval(hf_surf_var)) then
    write(0,*) trim(exec_name)//': total heat production ',hp_surf*hp_dep,' is greater than '//&
               'surface heat flow'
    call error_exit(1)
endif

! Propagate surface heat flow down to base and calculate temperature gradient
do j = 1,nt_total
    bas_grad(j) = (hf_surf_var(j)-hp_surf*hp_dep)*dz/cond_base
enddo




!**************************************************************************************************!
! *** Bulk Crustal Thickening/Thinning ************************************************************!
!**************************************************************************************************!
do i = 1,nthicken

    ! Check whether we are thickening or thinning
    if (thicken_dat(i,3).lt.0) then
        isThinning = .true.
    else
        isThinning = .false.
    endif
    thicken_dat(i,3) = abs(thicken_dat(i,3))

    ! Make sure thinning amount is less than crustal thickness
    if (isThinning .and. thicken_dat(i,3).ge.thicken_dat(i,5)) then
        write(0,*) trim(exec_name)//': thinning magnitude must be less than crustal thickness'
        write(0,*) 'thinning amount   =',thicken_dat(i,3)
        write(0,*) 'crustal thickness =',thicken_dat(i,5)
        call error_exit(1)
    endif

    ! Set thickening/thinning timesteps
    nstart = int(thicken_dat(i,1)/dt)             ! Starting timestep                          NN(1)
    nduration = int(thicken_dat(i,2)/dt)          ! Duration in timesteps                      NN(2)
    nthick = int(thicken_dat(i,3)/dz)             ! Amount to thicken/thin crust in nodes      NN(4)
    rate = dble(nthick)/dble(nduration)           ! Rate in nodes/timestep
    jbeg = nstart + 1                             ! First timestep of thickening/thinning
    jend = nstart + nduration                     ! Last timestep of thickening/thinning
    ct = 0                                        ! Initialize counter for thickening/thinning steps
    do j = jbeg,jend
        arg = dble(j-nstart)*rate - dble(ct)      ! Test for thickening/thinning at this timestep
        if (arg.ge.1.0d0) then

            if (isThinning) then
                action_array(j) = action_array(j) + 16
            else
                action_array(j) = action_array(j) + 8
            endif
            ct = ct + 1

            ! Save some parameters at beginning of event for bookkeeping
            if (thicken_start(i).eq.0) then
                thicken_start(i) = j                  ! Timestep to start event
                itop = int(thicken_dat(i,4)/dz) + 1   ! Top node of crust
                ithick = int(thicken_dat(i,5)/dz)     ! Initial crustal thickness in nodes
                crust_dat(i,1) = itop                 ! Save top node of crust
                crust_dat(i,2) = itop + ithick        ! Save initial bottom node of crust
            endif
        endif
    enddo
    if (isThinning) then
        thicken_dat(i,3) = -thicken_dat(i,3)
    endif

    ! Check whether tracked horizons should move during thickening/thinning
    if (thickenHorizons) then
        do ihor = 1,nhorizons
            ndep = depth_node(ihor)   ! horizon node

            ! How fast each horizon changes depth
            if (crust_dat(i,1).le.ndep .and. ndep.le.crust_dat(i,2)) then
                thickening_ratio = dble(nthick)*dz/thicken_dat(i,5)
                delta_depth = dble(ndep-itop-total_shift(ihor))*thickening_ratio
                rate = delta_depth/dble(nduration)
            elseif (crust_dat(i,2).lt.ndep) then
                rate = dble(nthick)/dble(nduration)
            else
                rate = 0.0d0
            endif

            ct = 0
            do j = jbeg,jend
                arg = dble(j-nstart)*rate - dble(ct)  ! Test for horizon moving

                if (arg.gt.1.0d0) then
                    if (isThinning) then
                        horizon_shift(ihor,j) = -1    ! Move horizon up by one node
                        total_shift(ihor) = total_shift(ihor) - 1
                    else
                        horizon_shift(ihor,j) = 1     ! Move horizon down by one node
                        total_shift(ihor) = total_shift(ihor) + 1
                    endif
                    ct = ct + 1
                endif
            enddo
        enddo
    endif
enddo




!**************************************************************************************************!
! *** Surface Temperature Variations **************************************************************!
!**************************************************************************************************!

! Add step temperature changes to surface temperature
do i = 1,ntempsteps
    jbeg = int(temp_step_dat(i,1)/dt)
    jend = nt_total
    do j = jbeg,jend
        temp_surf_var(j) = temp_surf_var(j) + temp_step_dat(i,2)
    enddo
enddo

! Add ramp temperature changes to surface temperature
do i = 1,ntempramps
    jbeg = int(temp_ramp_dat(i,1)/dt)
    jend = jbeg + int(temp_ramp_dat(i,2)/dt)
    if (jend.gt.nt_total) then
        write(0,*) trim(exec_name)//': end of temperature ramp period goes beyond model duration'
        call error_exit(1)
    endif
    rate = temp_ramp_dat(i,3)/temp_ramp_dat(i,2)
    do j = jbeg,jend
        temp_surf_var(j) = temp_surf_var(j) + (j-jbeg)*dt*rate
    enddo
    do j = jend+1,nt_total
        temp_surf_var(j) = temp_surf_var(j) + temp_ramp_dat(i,3)
    enddo
enddo

! Add sinusoidal temperature variations to surface temperature
do i = 1,ntempcycles
    jbeg = int(temp_cycle_dat(i,1)/dt) + 1
    jend = jbeg + int(temp_cycle_dat(i,2)/dt) - 1
    amp = temp_cycle_dat(i,3)
    freq = temp_cycle_dat(i,4)
    if (jend.gt.nt_total) then
        write(0,*) jbeg,jend,nt_total
        write(0,*) trim(exec_name)//': end of temperature sinusoid period goes beyond model duration'
        call error_exit(1)
    endif
    do j = jbeg,jend
        temp_surf_var(j) = temp_surf_var(j) + amp*sin(2.0d0*pi*freq*dble(j-jbeg)*dt)
    enddo
enddo




! Print timing of actions to file
if (timing_file.ne.'') then
    open(unit=13,file=timing_file,status='unknown')
    do i = 1,nburial
        write(13,*) 'burial',i,burial_dat(i,1),burial_dat(i,1)+burial_dat(i,2)
    enddo
    do i = 1,nuplift
        write(13,*) 'uplift',i,uplift_dat(i,1),uplift_dat(i,1)+uplift_dat(i,2)
    enddo
    do i = 1,nthrust
        write(13,*) 'thrust',i,thrust_dat(i,1)
    enddo
    do i = 1,nthicken
        write(13,*) 'thicken',i,thicken_dat(i,1),thicken_dat(i,1)+thicken_dat(i,2)
    enddo
    do i = 1,nhfvars
        write(13,*) 'hfvar',i,hfvar(i,1)
    enddo
    do i = 1,ntempsteps
        write(13,*) 'temp_step',i,temp_step_dat(i,1)
    enddo
    do i = 1,ntempramps
        write(13,*) 'temp_ramp',i,temp_ramp_dat(i,1),temp_ramp_dat(i,1)+temp_ramp_dat(i,2)
    enddo
    do i = 1,ntempcycles
        write(13,*) 'temp_step',i,temp_cycle_dat(i,1),temp_cycle_dat(i,1)+temp_cycle_dat(i,2)
    enddo
    close(13)
endif



if (verbosity.ge.2) then
    write(*,*) 'setup_action_arrays: finished'
endif


return
end subroutine






!--------------------------------------------------------------------------------------------------!






subroutine init_action_arrays()
!----
! Allocate memory to tectonic action arrays and initialize their values
!----

use tqtec


implicit none


! Free memory from tectonic action arrays - these should not be allocated, but ¯\_(ツ)_/¯
if (allocated(action_array)) then
    deallocate(action_array)
endif
if (allocated(burial_cond)) then
    deallocate(burial_cond)
endif
if (allocated(bas_grad)) then
    deallocate(bas_grad)
endif
if (allocated(thrust_step)) then
    deallocate(thrust_step)
endif
if (allocated(crust_dat)) then
    deallocate(crust_dat)
endif
if (allocated(thicken_start)) then
    deallocate(thicken_start)
endif
if (allocated(horizon_shift)) then
    deallocate(horizon_shift)
endif


! Allocate memory to tectonic action arrays
! DO NOT NEED TO STORE CONDUCTIVITY AND BASE TEMP GRADIENT FOR ALL TIMESTEPS.
! THOSE ARRAYS ARE MOSTLY ZEROS!
allocate(temp_surf_var(nt_total))
allocate(action_array(nt_total))    ! action code at each timestep: 0=no action, >0=action
if (nburial.gt.0) then
    allocate(burial_cond(nt_total)) ! conductivity of material added to model during burial step
endif
allocate(bas_grad(nt_total))        ! temperature gradient at the base of the model
if (nthrust.gt.0) then
    allocate(thrust_step(nthrust))  ! timestep for each thrust action
endif
if (nthicken.gt.0) then
    allocate(crust_dat(nthicken,2))
    allocate(thicken_start(nthicken))
endif
if (nthicken.gt.0.and.thickenHorizons) then
    allocate(horizon_shift(nhorizons,nt_total))
endif


! Initialize arrays
action_array = 0
temp_surf_var = temp_surf
if (nburial.gt.0) then
    burial_cond = 0.0d0
endif
bas_grad = 0.0d0
if (nthrust.gt.0) then
    thrust_step = 0
endif
if (nthicken.gt.0) then
    crust_dat = 0
    thicken_start = 0
endif
if (nthicken.gt.0.and.thickenHorizons) then
    horizon_shift = 0
endif


end subroutine




!--------------------------------------------------------------------------------------------------!




subroutine check_actions()

use tqtec

implicit none


! Local variables
integer :: i
double precision :: thick_init
double precision :: thrust_dep
double precision :: thick_end


! Check thrust geometry is compatible with model setup
do i = 1,nthrust

    thick_init = int(thrust_dat(i,3)/dz)  ! Thickness prior to thrusting, in nodes
    thrust_dep = int(thrust_dat(i,4)/dz)  ! Depth of emplacement, in nodes
    thick_end = int(thrust_dat(i,5)/dz)   ! Final thickness of thrust sheet, in nodes

    if (thrust_dat(i,5).le.0.0d0) then
        write(0,*) 'Final thrust sheet thickness must be greater than 0'
        write(0,*) 'Final thickness:   ',thrust_dat(i,5)
        call error_exit(1)
    endif
    if (thick_init.lt.thick_end) then
        write(0,*) 'Final thrust sheet thickness must be <= initial thickness'
        write(0,*) 'Initial thickness: ',thrust_dat(i,3)
        write(0,*) 'Final thickness:   ',thrust_dat(i,5)
        call error_exit(1)
    endif
    if (thrust_dat(i,4).gt.thrust_dat(i,4)) then
        write(0,*) 'Thrust sheet emplacement depth must be less than or equal to final thickness'
        write(0,*) 'Emplacement depth: ',thrust_dat(i,4)
        write(0,*) 'Final thickness:   ',thrust_dat(i,5)
        call error_exit(1)
    endif
    if (2*thick_init.ge.nnodes) then
        write(0,*) 'tqtec: not enough nodes in model to smooth thrust sheet geotherm'
        write(0,*) 'thrust sheet is',thick_init,' nodes'
        write(0,*) 'model is       ',nnodes    ,' nodes'
        write(0,*) 'model must be > 2*thrust'
        call error_exit(1)
    endif

enddo



return
end subroutine








!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!------------------------------- TECTONIC ACTION ROUTINES -----------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!



subroutine bury()
!----
! Bury the horizons by shifting physical parameters down node list and updating the surface node
!----

use tqtec, only: verbosity, &
                 nnodes, &
                 istep, &
                 conductivity, &
                 temp, &
                 hp, &
                 nhorizons, &
                 depth_node, &
                 temp_surf, &
                 burial_cond

implicit none

! Local variables
integer :: i, j


if (verbosity.ge.3) then
    write(*,*) '    burying...'
endif


! Shift physical parameters (temperature, heat production, conductivity) down by one node
do i = 1,nnodes-1
    j = nnodes-i
    temp(j+1) = temp(j)
    hp(j+1) = hp(j)
    conductivity(j+1) = conductivity(j)
enddo

! Update the top node temperature, heat production, and conductivity
temp(1) = temp_surf
hp(1) = 0.0d0
conductivity(1) = burial_cond(istep)

! Move all of the tracked horizons down by one node
do i = 1,nhorizons
    depth_node(i) = depth_node(i)+1
enddo

return
end subroutine



!--------------------------------------------------------------------------------------------------!



subroutine erode()
!----
! Erode and uplift the horizons by shifting physical parameters up node list by one and removing
! the surface node
!----

use tqtec, only: verbosity, &
                 nnodes, &
                 temp, &
                 hp, &
                 conductivity, &
                 nhorizons, &
                 depth_node, &
                 cond_base, &
                 dtemp_wo_hp

implicit none

! Local variables
integer :: i


if (verbosity.ge.3) then
    write(*,*) '    eroding...'
endif


! Shift physical parameters (temperature, heat production, conductivity) up by one node
do i = 2,nnodes
    temp(i-1) = temp(i)
    hp(i-1) = hp(i)
    conductivity(i-1) = conductivity(i)
enddo

! Update the bottom node temperature, heat production, and conductivity
temp(nnodes) = temp(nnodes-1) + dtemp_wo_hp
hp(nnodes) = 0.0d0
conductivity(nnodes) = cond_base

! Move horizons upward by one node
do i = 1,nhorizons
    depth_node(i) = depth_node(i)-1
    ! if (depth_node(i).le.0) then
    !     depth_node(i) = 0
    ! endif
enddo

return
end subroutine



!--------------------------------------------------------------------------------------------------!



subroutine thrust(plate)
!----
! Generate a thrust fault, moving the horizons to the lower plate
!----

use tqtec, only: verbosity, &
                 nnodes, &
                 istep, &
                 dz, &
                 conductivity, &
                 temp, &
                 hp, &
                 nhorizons, &
                 depth_node, &
                 temp_surf, &
                 nthrust, &
                 thrust_step, &
                 thrust_dat


implicit none


! Arguments
character(len=*) :: plate


! Local variables
integer :: i, k, ismooth, nsmooth
integer :: thick_init, thrust_dep, thick_end, ierosion, dnode
double precision :: upl_conductivity(nnodes), upl_hp(nnodes), upl_temp(nnodes)



if (verbosity.ge.3) then
    write(*,*) '    thrusting horizons into ',trim(plate),' plate...'
endif


! Set the number of smoothing passes
nsmooth = 10


! Set the thrust number
k = 0
do i = 1,nthrust
    if (thrust_step(i).eq.istep) then
        k = i
        exit
    endif
enddo
if (k.le.0.or.k.gt.nthrust) then
    write(0,*) 'Error setting thrust number at istep=',istep
    call error_exit(1)
endif


! Save thrust sheet geometric parameters
thick_init = int(thrust_dat(k,3)/dz)  ! Thickness of thrust sheet prior to thrusting, in nodes
thrust_dep = int(thrust_dat(k,4)/dz)  ! Depth of emplacement, in nodes
thick_end = int(thrust_dat(k,5)/dz)   ! Final thickness of thrust sheet, in nodes


! Number of nodes to be instantaneously removed from the top of the thrust sheet
ierosion = thick_init - thick_end

! Save the properties of nodes in the thrust sheet, minus the nodes eroded off the top
do i = ierosion+1,thick_init
    upl_conductivity(i-ierosion) = conductivity(i)
    upl_hp(i-ierosion) = hp(i)
    upl_temp(i-ierosion) = temp(i)
enddo


! Copy the existing model array (the lower plate/footwall) from the depth of thrust sheet
! emplacement to the base of the model to the top of the model. This removes the array 
! representing the lower plate/footwall from the surface down to the depth of emplacement.
do i = thrust_dep+1,nnodes
    conductivity(i-thrust_dep) = conductivity(i)
    hp(i-thrust_dep) = hp(i)
    temp(i-thrust_dep) = temp(i)
enddo


! Shift the model array (lower plate/footwall) down by the final thickness of the thrust sheet
! to make room for the duplicated nodes (upper plate/hanging wall)
do i = nnodes,thick_end+1,-1
    conductivity(i) = conductivity(i-thick_end)
    hp(i) = hp(i-thick_end)
    temp(i) = temp(i-thick_end)
enddo


! Put the thrust sheet nodes (upper plate/hanging wall) on top
do i = 1,thick_end
    conductivity(i) = upl_conductivity(i)
    hp(i) = upl_hp(i)
    temp(i) = upl_temp(i)
enddo


! Move tracked horizons into correct location in upper or lower plate
if (plate.eq.'upper') then
    do i = 1,nhorizons
        depth_node(i) = depth_node(i) - ierosion
        if (depth_node(i).le.0) then
            depth_node(i) = 0
        endif
    enddo
elseif (plate.eq.'lower') then
    dnode = thick_end - thrust_dep
    do i = 1,nhorizons
        depth_node(i) = depth_node(i) + dnode
    enddo
else
    write(0,*) 'no plate named "',trim(plate),'"'
    call error_exit(1)
endif


! Smooth temperatures where there is a sharp thermal gradient due to instantaneous thrusting
do ismooth = 1,nsmooth
    temp(1) = (temp_surf+temp(2))/2.0d0
    temp(2) = (temp(1)+temp(2)+temp(3))/3.0d0
    do i = 3,2*thick_init
        temp(i) = (temp(i-2)+temp(i-1)+temp(i)+temp(i+1)+temp(i+2))/5.0d0
    enddo
enddo

return
end subroutine



!--------------------------------------------------------------------------------------------------!



subroutine thicken()

use tqtec, only: nnodes, &
                 conductivity, &
                 temp, &
                 crust_top, &
                 crust_bot, &
                 verbosity

implicit none

! Local variables
integer :: i
integer :: j
integer :: crust_thick
double precision :: ratio


if (verbosity.ge.3) then
    write(*,*) '    thickening...'
endif


! Duplicate conductivity at new crust node and shift conductivities down one node
conductivity(crust_bot+1:nnodes) = conductivity(crust_bot:nnodes-1)

! Duplicate temperature at new crust node and shift temperatures down one node
temp(crust_bot+1:nnodes) = temp(crust_bot:nnodes-1)


! Redistribute temperatures throughout thickened crust
crust_thick = crust_bot - crust_top
ratio = dble(crust_thick)/dble(crust_thick+1)
do i = 0,crust_thick
    j = i + crust_top
    temp(j) = temp(j) * (1 + dble(i)*(ratio-1.0d0)/dble(crust_thick+1))
enddo


! Count the duplicated node at the bottom of the crust
crust_bot = crust_bot + 1


return
end subroutine


!--------------------------------------------------------------------------------------------------!



subroutine thin()

use tqtec, only: nnodes, &
                 conductivity, &
                 temp, &
                 crust_top, &
                 crust_bot, &
                 verbosity

implicit none

! Local variables
integer :: i
integer :: j
integer :: crust_thick
double precision :: ratio


if (verbosity.ge.3) then
    write(*,*) '    thinning...'
endif


! Remove bottom crust node and shift conductivities up one node
conductivity(crust_bot:nnodes-1) = conductivity(crust_bot+1:nnodes)

! Shift temperatures up one node
temp(crust_bot:nnodes-1) = temp(crust_bot+1:nnodes)


! Redistribute temperatures throughout thinned crust
crust_thick = crust_bot - crust_top
ratio = dble(crust_thick)/dble(crust_thick-1)
do i = 0,crust_thick
    j = i + crust_top
    temp(j) = temp(j) * (1 + dble(i)*(ratio-1.0d0)/dble(crust_thick-1))
enddo


! Remove a node from the bottom of the crust
crust_bot = crust_bot - 1


return
end subroutine


