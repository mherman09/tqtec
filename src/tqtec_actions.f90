!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!----------------------------------- TECTONIC ACTIONS ---------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!


subroutine setup_action_arrays()
!----
! Define arrays to control tectonic action events
!----

use tqtec, only: timing_file, &
                 verbosity, &
                 nt_total, &
                 dz, &
                 dt, &
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
                 thicken_dat

implicit none


! Local variables
integer :: i, ct
integer :: j, jbeg, jend
integer :: nstart, nduration, nthick, top_crust, crust_thick, bot_crust
double precision :: rate, arg
double precision :: hf_surf_var(nt_total)
integer :: intqt(nhfvars)



if (verbosity.ge.2) then
    write(*,*) 'setup_action_arrays: starting'
endif


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


! Allocate memory to tectonic action arrays
allocate(action(nt_total))         ! action code at each timestep: 0=no action, >0=action
allocate(burial_cond(nt_total))    ! conductivity of material added to model during burial step
allocate(bas_grad(nt_total))       ! temperature gradient at the base of the model
if (nthrust.gt.0) then
    allocate(thrust_step(nthrust)) ! thrust timestep for each thrust action
endif


! Initialize arrays
action = 0
burial_cond = 0.0d0
bas_grad = 0.0d0
if (nthrust.gt.0) then
    thrust_step = 0
endif


! Burial periods
do i = 1,nburial
    nstart = int(burial_dat(i,1)/dt)              ! Starting timestep
    nduration = int(burial_dat(i,2)/dt)           ! Duration in timesteps
    nthick = int(burial_dat(i,3)/dz)              ! Thickness in nodes
    rate = dble(nthick)/dble(nduration)           ! Rate in nodes/timestep
    jbeg = nstart + 1                             ! First timestep
    jend = nstart + nduration                     ! Last timestep
    ct = 0                                        ! Initialize counter for number of burial increments
    do j = jbeg,jend
        arg = dble(j-nstart)*rate - dble(ct)      ! Test for burying at this timestep
        if (arg.ge.1.0d0) then
            action(j) = 1                         ! Set burial action code = 1
            burial_cond(j) = burial_dat(i,4)
            ct = ct + 1
        endif
    enddo
enddo


! Uplift/erosion periods
do i = 1,nuplift
    nstart = int(uplift_dat(i,1)/dt)              ! Starting timestep
    nduration = int(uplift_dat(i,2)/dt)           ! Duration in timesteps
    nthick = int(uplift_dat(i,3)/dz)              ! Thickness in nodes
    rate = dble(nthick)/dble(nduration)           ! Rate in nodes/timestep
    jbeg = nstart + 1                             ! First timestep
    jend = nstart + nduration                     ! Last timestep
    ct = 0                                        ! Initialize counter for number of uplift increments
    do j = jbeg,jend
        arg = dble(j-nstart)*rate - dble(ct)      ! Test for eroding at this timestep
        if (arg.ge.1.0d0) then
            action(j) = 2                         ! Set uplift action code = 2
            ct = ct + 1
        endif
    enddo
enddo


! Thrust events
do i = 1,nthrust
    nstart = int(thrust_dat(i,1)/dt) + 1          ! Timestep of thrust faulting
    action(nstart) = 3                            ! Set thrust action code = 3
    thrust_step(i) = nstart                       ! Save thrust timestep for each thrust action
    ! THTYPE(I): thrust_dat(i,2)
enddo


! Basal heat flow
! Initialize heat flow over time to be surface heat flow
hf_surf_var = hf_surf

! Timing of heat flow changes
do i = 1,nhfvars
    intqt(i) = int(hfvar(i,1)/dt)                 ! Timestep of heat flow change
enddo
do i = 1,nhfvars-1
    jbeg = intqt(i)
    jend = intqt(i+1)
    do j = jbeg,jend
        hf_surf_var(j) = hfvar(i,2)
    enddo
enddo
if (nhfvars.ge.1) then
    jbeg = intqt(nhfvars)
    jend = nt_total
    do j = jbeg,jend
        hf_surf_var(j) = hfvar(nhfvars,2)
    enddo
endif

! Check that heat production is never greater than surface heat flow
if (hp_surf*hp_dep.ge.maxval(hf_surf_var)) then
    write(0,*) 'tqtec: total heat production ',hp_surf*hp_dep,' is greater than surface heat flow'
    stop 1
endif

! Propagate surface heat flow down to base and calculate temperature gradient
do j = 1,nt_total
    bas_grad(j) = (hf_surf_var(j)-hp_surf*hp_dep)*dz/cond_base
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
    close(13)
endif


if (verbosity.ge.2) then
    write(*,*) 'setup_action_arrays: finished'
endif

return
end subroutine


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