!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------- PREPARE ACTION ARRAYS ------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!


subroutine setup_action_arrays()
!----
! Construct arrays to control tectonic action events
!----

use tqtec, only: timing_file, &
                 verbosity, &
                 nt_total, &
                 dz, &
                 dt, &
                 nhorizons, &
                 depth_node, &
                 cond_base, &
                 hf_surf, &
                 hp_surf, &
                 hp_dep, &
                 nburial, &
                 burial_dat, &
                 nuplift, &
                 uplift_dat, &
                 nthrust, &
                 thrust_dat, &
                 thrust_step, &
                 nhfvars, &
                 hfvar, &
                 bas_grad, &
                 action, &
                 burial_cond, &
                 nthicken, &
                 thicken_dat, &
                 thicken_start, &
                 crust_dat, &
                 thickenHorizons, &
                 horizon_shift

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
logical :: isThinning



if (verbosity.ge.2) then
    write(*,*) 'setup_action_arrays: starting'
endif



! Initialize global tectonic action arrays
call init_action_arrays()



! Initialize local tectonic action arrays
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
! *** Burial Events ***
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
            action(j) = 1                         ! Set burial action code = 1
            burial_cond(j) = burial_dat(i,4)
            ct = ct + 1
        endif
    enddo
enddo



!**************************************************************************************************!
! *** Uplift/Erosion Events ***********************************************************************!
do i = 1,nuplift
    nstart = int(uplift_dat(i,1)/dt)              ! Starting timestep
    nduration = int(uplift_dat(i,2)/dt)           ! Duration in timesteps
    nthick = int(uplift_dat(i,3)/dz)              ! Thickness in nodes
    rate = dble(nthick)/dble(nduration)           ! Rate in nodes/timestep
    jbeg = nstart + 1                             ! First timestep
    jend = nstart + nduration                     ! Last timestep
    ct = 0                                        ! Initialize counter for uplift increments
    do j = jbeg,jend
        arg = dble(j-nstart)*rate - dble(ct)      ! Test for eroding at this timestep
        if (arg.ge.1.0d0) then
            action(j) = 2                         ! Set uplift action code = 2
            ct = ct + 1
        endif
    enddo
enddo



!**************************************************************************************************!
! *** Thrust Events *******************************************************************************!
do i = 1,nthrust
    nstart = int(thrust_dat(i,1)/dt) + 1          ! Timestep of thrust faulting
    action(nstart) = 3                            ! Set thrust action code = 3
    thrust_step(i) = nstart                       ! Save thrust timestep for each thrust action
    ! THTYPE(I): thrust_dat(i,2)
enddo


!**************************************************************************************************!
! *** Surface Heat Flow Changes *******************************************************************!
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
    write(0,*) 'tqtec: total heat production ',hp_surf*hp_dep,' is greater than surface heat flow'
    call error_exit(1)
endif

! Propagate surface heat flow down to base and calculate temperature gradient
do j = 1,nt_total
    bas_grad(j) = (hf_surf_var(j)-hp_surf*hp_dep)*dz/cond_base
enddo




!**************************************************************************************************!
! *** Bulk Crustal Thickening/Thinning ************************************************************!
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
        write(0,*) 'tqtec: thinning magnitude must be less than crustal thickness'
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
                action(j) = -4                    ! Set thinning action code = -4
            else
                action(j) = 4                     ! Set thickening action code = 4
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



! Print timing of tectonic actions to file
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

use tqtec, only: nt_total, &
                 nhorizons, &
                 action, &
                 nburial, &
                 burial_cond, &
                 bas_grad, &
                 nthrust, &
                 thrust_step, &
                 crust_dat, &
                 nthicken, &
                 thicken_start, &
                 thickenHorizons, &
                 horizon_shift


implicit none


! Free memory from tectonic action arrays - these should not be allocated, but ¯\_(ツ)_/¯
if (allocated(action)) then
    deallocate(action)
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
allocate(action(nt_total))         ! action code at each timestep: 0=no action, >0=action
if (nburial.gt.0) then
    allocate(burial_cond(nt_total))! conductivity of material added to model during burial step
endif
allocate(bas_grad(nt_total))       ! temperature gradient at the base of the model
if (nthrust.gt.0) then
    allocate(thrust_step(nthrust)) ! timestep for each thrust action
endif
if (nthicken.gt.0) then
    allocate(crust_dat(nthicken,2))
    allocate(thicken_start(nthicken))
endif
if (nthicken.gt.0.and.thickenHorizons) then
    allocate(horizon_shift(nhorizons,nt_total))
endif


! Initialize arrays
action = 0
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

    if (thick_init.lt.thick_end) then
        write(0,*) 'tqtec: final thrust sheet thickness must be <= initial thickness'
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



subroutine thrust_upperplate()
!----
! Generate a thrust fault, keeping the horizons in the upper plate
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

! Local variables
integer :: i, k, ismooth, nsmooth
integer :: thick_init, thrust_dep, thick_end, ierosion
double precision :: upl_conductivity(nnodes), upl_hp(nnodes), upl_temp(nnodes)


if (verbosity.ge.3) then
    write(*,*) '    thrusting horizons into upper plate...'
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

! Save thrust sheet parameters
thick_init = int(thrust_dat(k,3)/dz)  ! Thickness prior to thrusting, in nodes
thrust_dep = int(thrust_dat(k,4)/dz)  ! Depth of emplacement, in nodes
thick_end = int(thrust_dat(k,5)/dz)   ! Final thickness of thrust sheet, in nodes


! C     COPY THE PART OF THE ARRAY THAT WILL BE THRUSTED AND ERODE OFF
! C     AMOUNT THAT GETS ERODED DURING THRUST EVENT.
ierosion = thick_init - thick_end
do i = ierosion+1,thick_init
    upl_conductivity(i-ierosion) = conductivity(i)
    upl_hp(i-ierosion) = hp(i)
    upl_temp(i-ierosion) = temp(i)
enddo

! C     REMOVE THE OLD ARRAY (LOWER PLATE) DOWN TO THE DEPTH OF
! C     EMPLACEMENT AND MOVE THE REST OF THE ARRAY DOWN TO MAKE ROOM
! C     FOR THE UPPER PLATE.
do i = thrust_dep+1,nnodes
    conductivity(i-thrust_dep) = conductivity(i)
    hp(i-thrust_dep) = hp(i)
    temp(i-thrust_dep) = temp(i)
enddo
do i = nnodes,thick_end+1,-1
    conductivity(i) = conductivity(i-thick_end)
    hp(i) = hp(i-thick_end)
    temp(i) = temp(i-thick_end)
enddo

! C     PUT THE TWO ARRAYS TOGETHER
! I.e., put the thrust sheet on top
do i = 1,thick_end
    conductivity(i) = upl_conductivity(i)
    hp(i) = upl_hp(i)
    temp(i) = upl_temp(i)
enddo

! C     MOVE POINTS OF INTEREST AROUND FOR UPPER PLATE:
! I.e., move the specified horizons into the upper plate
do i = 1,nhorizons
    depth_node(i) = depth_node(i) - ierosion
    if (depth_node(i).le.0) then
        depth_node(i) = 0
    endif
enddo
! WHAT HAPPENS IF THE THRUST SHEET IS THINNER THAN THE DEEPEST HORIZON? THIS HORIZON CANNOT GO INTO
! THE UPPER PLATE THEN...

! Smooth temperatures where there is a sharp thermal gradient due to instantaneous thrusting
do ismooth = 1,nsmooth
    temp(1) = (temp_surf+temp(2))/2.0d0
    temp(2) = (temp(1)+temp(2)+temp(3))/3.0d0
    do i = 3,2*thick_init
        temp(i) = (temp(i-2)+temp(i-1)+temp(i)+temp(i+1)+temp(i+2))/5.0d0
    enddo
enddo
! NOTE: Smoothing does not take into account a thrust sheet emplaced at depth

return
end subroutine



!--------------------------------------------------------------------------------------------------!



subroutine thrust_lowerplate()
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

! Local variables
integer :: i, k, ismooth, nsmooth
integer :: thick_init, thrust_dep, thick_end, ierosion, dnode
double precision :: upl_conductivity(nnodes), upl_hp(nnodes), upl_temp(nnodes)


if (verbosity.ge.3) then
    write(*,*) '    thrusting horizons into lower plate...'
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

! Save thrust sheet parameters
thick_init = int(thrust_dat(k,3)/dz)  ! Thickness prior to thrusting, in nodes
thrust_dep = int(thrust_dat(k,4)/dz)  ! Depth of emplacement, in nodes
thick_end = int(thrust_dat(k,5)/dz)   ! Final thickness of thrust sheet, in nodes


dnode = thick_end - thrust_dep

! C     COPY THE PART OF THE ARRAY THAT WILL BE THRUSTED AND ERODE OFF
! C     THE AMOUNT THAT GETS ERODED DURING THE THRUST EVENT.
ierosion = thick_init - thick_end
do i = ierosion+1,thick_init
    upl_conductivity(i-ierosion) = conductivity(i)
    upl_hp(i-ierosion) = hp(i)
    upl_temp(i-ierosion) = temp(i)
enddo

! C     REMOVE THE OLD ARRAY (LOWER PLATE) DOWN TO THE DEPTH OF
! C     EMPLACEMENT AND MOVE THE REST OF THE ARRAY DOWN TO MAKE ROOM
! C     FOR THE UPPER PLATE.
do i = thrust_dep+1,nnodes
    conductivity(i-thrust_dep) = conductivity(i)
    hp(i-thrust_dep) = hp(i)
    temp(i-thrust_dep) = temp(i)
enddo
do i = nnodes,thick_end+1,-1
    conductivity(i) = conductivity(i-thick_end)
    hp(i) = hp(i-thick_end)
    temp(i) = temp(i-thick_end)
enddo

! C     PUT THE TWO ARRAYS TOGETHER
! I.e., put the thrust sheet on top
do i = 1,thick_end
    conductivity(i) = upl_conductivity(i)
    hp(i) = upl_hp(i)
    temp(i) = upl_temp(i)
enddo

! C     MOVE POINTS OF INTEREST AROUND
! C     FOR LOWER PLATE:
do i = 1,nhorizons
    depth_node(i) = depth_node(i) + dnode
enddo

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


