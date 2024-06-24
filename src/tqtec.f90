!----
! TQTec (Temperature, Heat Flow, Tectonics)
!
! Authors:
!     - Kevin Furlong (original Fortran 77 program)
!     - Matt Herman (Modern Fortran version, i.e., what you are looking at right now!)
!
! C     CALCULATES THE ONE-DIMENSIONAL TRANSIENT THERMAL FIELD WITHIN
! C     AN AREA THAT UNDERGOES EPISODES OF:  BURIAL, EROSION, AND
! C     THRUSTING.  THIS PROGRAM CAN BE USED TO MONITOR POINTS THAT
! C     MOVE FROM THE UPPER TO LOWER (or v.v) PLATES.
!
! Incorporates (or will eventually incorporate) bulk thickening and thinning capabilities
! done by Chris Guzofski as part of his Master's thesis
!----


module tqtec


! Inputs/outputs
character(len=512) :: input_file                  ! name of input file                              INFIL
character(len=8) :: input_mode                    ! how to read input parameters (user, file)
character(len=512) :: output_file                 ! name of output file                             OUTFIL
character(len=512) :: geotherm_file               ! name of output geotherm file
character(len=512) :: timing_file                 ! name of output tectonic action timing file
integer :: verbosity                              ! level of information to print during execution


! Finite difference parameters
integer :: nnodes                                 ! number of spatial nodes                         N
integer :: nt_total                               ! number of time steps                            Q1 (updated), II(5)
integer :: istep                                  ! current time step                               V
double precision :: dz                            ! node spacing (km)                               H1, II(1)
double precision :: dt                            ! time step interval (Ma)                         K1, II(2)
double precision :: r1                            ! finite difference time factor                   R1


! Timing
double precision :: t_total                       ! total model time (Ma)                           Q1 (initial), II(5)
double precision :: t_geotherm_output             ! time per geotherm output (Ma)                   M1 (initial)
integer :: nt_geotherm_output                     ! time steps between geotherm outputs             M1 (updated)


! Nodal parameters
double precision, allocatable :: conductivity(:)  ! conductivity                                    COND
double precision, allocatable :: temp(:)          ! temperature                                     B
double precision, allocatable :: hp(:)            ! heat production                                 H
double precision, allocatable :: hf(:)            ! heat flow                                       Q


! Horizons of interest
integer :: nhorizons                              ! number of horizons                              10
double precision, allocatable :: depth(:)         ! depth of horizons                               Y (initial)
integer, allocatable :: depth_node(:)             ! horizon nodes                                   Y (updated)


! Material properties
integer :: nlayers                                ! number of distinct material layers              INL
double precision, allocatable :: layer(:,:)       ! layer(:,1): depth to top (km)                   TOP
                                                  ! layer(:,2): thickness (km)                      THICK
                                                  ! layer(:,3): conductivity (W/(m*K))              ACOND
double precision :: diffusivity                   ! diffusivity                                     D1, II(6)
double precision :: cond_base                     ! basal conductivity                              C1


! Boundary conditions
double precision :: temp_surf                     ! surface temperature (C)                         W(1)
double precision :: hp_surf                       ! surface heat production                         A1, II(3)
double precision :: hp_dep                        ! depth of heat production                        B1, II(4)
double precision :: hf_surf                       ! surface heat flow                               G1
double precision :: hf_base                       ! basal heat flow                                 QBASE
double precision :: dtemp_wo_hp                   ! temp change w/o heat prod                       W(2)
double precision :: temp_factor                   ! temp scaling factor                             W1
double precision :: temp_base_adj                 ! temp at node nnodes+1                           W(3)


! Tectonic events
integer, allocatable :: action(:)                 ! burial (1), erosion (2), thrust (3)             P
logical :: isTectonicActionDefined

integer :: nburial !------------------------------! number of burial events                         NBP
double precision, allocatable :: burial_dat(:,:)  ! burial_dat(:,1): start (Ma)                     AN(1)
                                                  ! burial_dat(:,2): duration (Ma)                  AN(2)
                                                  ! burial_dat(:,3): thickness (km)                 AN(3)
                                                  ! burial_dat(:,4): conductivity (W/(m*K))         AN(4)
double precision, allocatable :: burial_cond(:)   ! conductivity of new material at each timestep   BCOND

integer :: nuplift !------------------------------! number of uplift/erosion events                 NUEP
double precision, allocatable :: uplift_dat(:,:)  ! uplift_dat(:,1): start (Ma)                     AN(1)
                                                  ! uplift_dat(:,2): duration (Ma)                  AN(2)
                                                  ! uplift_dat(:,3): thickness (km)                 AN(3)

integer :: nthrust !------------------------------! number of thrusting events                      NTP
integer, allocatable :: thrust_step(:)            ! thrusting timestep array
double precision, allocatable :: thrust_dat(:,:)  ! thrust_dat(:,1): start (Ma)                     AN(1)
                                                  ! thrust_dat(:,2): upper (1) or lower (2) plate   AN(2)
                                                  ! thrust_dat(:,3): initial base (km)              AZ(1)
                                                  ! thrust_dat(:,4): initial depth (km)             AZ(2)
                                                  ! thrust_dat(:,5): initial thickness (km)         AZ(3)

integer :: nhfvars !------------------------------! number of surface heat flow variations          QSTEP
double precision, allocatable :: hfvar(:,:)       ! (1) start (2) new heat flow                     QVTIME
double precision, allocatable :: bas_grad(:)      ! temperature gradient at the base of the model   BASGRAD

integer :: nthicken !-----------------------------! number of thickening/thinning events
double precision, allocatable :: thicken_dat(:,:) ! thicken_dat(:,1): start (Ma)                    AN(1)
                                                  ! thicken_dat(:,2): duration (Ma)                 AN(2)
                                                  ! thicken_dat(:,3): total thickening amount (km)  AN(4)
                                                  ! thicken_dat(:,4): top of crust (km)
                                                  ! thicken_dat(:,5): initial crustal thickness (km)
integer, allocatable :: thicken_start(:)          ! thicken_start(:): timestep thickening event starts
integer, allocatable :: crust_dat(:,:)            ! crust_dat(:,1): top of crust (node)
                                                  ! crust_dat(:,2): bottom of crust (node)
integer:: crust_top
integer:: crust_bot
logical :: thickenHorizons                        ! T: move tracked horizon depths with thickening


! Results array
double precision, allocatable :: results(:,:,:)   ! temperature and depth for each timstep/depth    r1


end module tqtec




!==================================================================================================!
!==================================================================================================!
!==================================================================================================!
!==================================================================================================!
!==================================================================================================!




program main

!----
! Solve for the 1-D transient thermal field defined by boundary conditions, material properties,
! and tectonic/geologic actions
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



! Initialize default model parameters
call initialize_defaults()                            ! tqtec_io.f90



! Parse command line options
call gcmdln()                                         ! tqtec_io.f90



if (verbosity.ge.1) then
    write(*,*) 'tqtec: starting'
endif



! Read control file or user input from standard input (formerly INPUT)
call read_model_parameters()                          ! tqtec_io.f90



! Check that model is deep enough to include all horizons of interest
if (dz*dble(nnodes).lt.maxval(depth)) then
    write(0,*) 'tqtec: model extent is shallower than deepest horizon'
    call error_exit(1)
endif



! Calculate model parameters
nt_total = int(t_total/dt)                            ! number of time steps
nt_geotherm_output = int(t_geotherm_output/dt)        ! number of time steps between geotherm outputs
temp_factor = diffusivity*dt/cond_base                ! temperature scaling factor
r1 = diffusivity*dt/(dz*dz)                           ! finite difference time factor
dtemp_wo_hp = (hf_surf - hp_surf*hp_dep)*dz/cond_base ! temperature difference without heat production



! Allocate arrays
allocate(depth_node(nhorizons))                       ! Node where each horizon of interest sits
allocate(conductivity(nnodes))                        ! Conductivity at each node
allocate(temp(nnodes))                                ! Temperature at each node
allocate(hp(nnodes))                                  ! Heat production at each node
allocate(hf(nt_total))                                ! Surface heat flow at each time
allocate(results(nt_total,2,nhorizons))               ! Temperature and depth at each time step for each horizon



! Set up tectonic action timing arrays (formerly HIST)
call setup_action_arrays()                            ! tqtec_actions.f90
call check_actions()                                  ! tqtec_actions.f90



! Initialize the temperature, heat flow, and heat production at each node (formerly INIT)
call initialize_thermal_parameters()                  ! tqtec.f90



! Print model parameters to standard output
if (verbosity.ge.1) then
    call print_model_parameters()                     ! tqtec_io.f90
endif



! Print the geotherm to a file (-geotherm flag)
if (geotherm_file.ne.'') then
    open(unit=12,file=geotherm_file,status='unknown')
    ! Header contains time step, time since start in Ma, and time until end in Ma
    write(12,'(A,I10,2F10.3)') '> #',0,0.0d0,0.0d0-t_total
    do i = 1,nnodes
        write(12,*) temp(i),dble(i)*dz
    enddo
endif



! Step through time and run the finite difference procedure to calculate the temperature at each
! node for each time step
istep = 0
do while (istep.lt.nt_total)


    ! Update the adjusted temperature at the base of the model
    temp_base_adj = temp(nnodes) + bas_grad(istep+1)


    ! Calculate the updated temperatures at each node (***THE MAIN FINITE DIFFERENCE PROCEDURE***)
    ! (formerly MAT and TRID)
    call update_temps(nnodes,ierr)                    ! tqtec.f90
    if (ierr.ne.0) then
        write(0,*) 'tqtec: error in update_temps() TRID algorithm at step',istep
        call error_exit(1)
    endif


    ! Increment the time step
    istep = istep + 1
    if (verbosity.le.2) then
        if (istep.lt.nt_total) then
            write(*,'(A,I6,A,I6,A)',advance='no') ' tqtec: working on step',istep,' of',nt_total,char(13)
        else
            write(*,'(A,I6,A,I6)') ' tqtec: working on step',istep,' of',nt_total
        endif
    elseif (verbosity.ge.3) then
        if (istep.lt.nt_total.and.action(istep).eq.0) then
            write(*,'(A,I6,A,I6,A)',advance='no') ' tqtec: working on step',istep,' of',nt_total,char(13)
        else
            write(*,'(A,I6,A,I6)') ' tqtec: working on step',istep,' of',nt_total
        endif
    endif


    ! Tectonic actions (in tqtec_actions.f90)
    if (action(istep).eq.1) then
        ! Add material to the top of the model (burial)
        call bury() ! (Formerly: BURIAL)

    elseif (action(istep).eq.2) then
        ! Remove material from the top of the model (erosion)
        call erode() ! (Formerly: EROS)

    elseif (action(istep).eq.3) then
        ! Get thrust event number
        j = 0
        do i = 1,nthrust
            if (thrust_step(i).eq.istep) then
                j = i
                exit
            endif
        enddo
        if (j.eq.0) then
            write(0,*) 'tqtec: error setting thrust number'
        endif
        iplate = int(thrust_dat(j,2))
        if (iplate.eq.0) then
            write(0,*) 'tqtec: error setting thrust number'
        endif
        ! Duplicate material and insert it (thrusting)
        if (iplate.eq.1) then
            ! Track horizons in the hanging wall (upper plate) of the thrust sheet
            call thrust_upperplate() ! (Formerly: THSTUP)
        elseif (iplate.eq.2) then
            ! Track horizons in the footwall (lower plate) of the thrust sheet
            call thrust_lowerplate() ! (Formerly: THSTLP)
        endif

    elseif (action(istep).eq.4) then
        ! Check to see if this is a new crustal thickening/thinning event
        do i = 1,nthicken
            ! write(0,*) 
            ! write(0,*) 'istep:              ',thicken_start(i)
            ! write(0,*) 'thicken_start:      ',thicken_start(i)
            if (thicken_start(i).eq.istep) then
                ! Reset initial crustal top and bottom nodes
                crust_top = crust_dat(i,1)
                crust_bot = crust_dat(i,2)
                write(0,*) 'crust_top:',crust_top
                write(0,*) 'crust_bot:',crust_bot
                exit
            endif
        enddo
        call thicken()
        write(*,*) 'temp(crust_bot):',temp(crust_bot)
        write(*,*) 'temp(nnodes):   ',temp(nnodes)

    elseif (action(istep).ne.0) then
        write(0,*) ''
        write(0,*) 'tqtec: error defining action flag'
        write(0,*) 'No action flag number',action(istep)
        call error_exit(1)

    endif


    ! Calculate surface heat flow for this time step
    hf(istep) = (temp(10)-temp(5))/(5.0d0*dz)   ! Temperature gradient near surface
    cond_surf = 0.0d0                           ! Average surface conductivity
    do i = 1,5
        cond_surf = cond_surf + conductivity(i+4)
    enddo
    cond_surf = cond_surf/5.0d0
    hf(istep) = hf(istep)*cond_surf             ! Heat flow = dT/dz * conductivity


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
            ! Header contains time step, time since start in Ma, and time until end in Ma
            write(12,'(A,I10,2F10.3)') '> #',istep,dble(istep)*dt,dble(istep)*dt-t_total
            do i = 1,nnodes
                write(12,*) temp(i),dble(i)*dz
            enddo
        endif
    endif

enddo   ! end of timestep loop


! Close the geotherm file if needed
if (geotherm_file.ne.'') then
    close(12)
endif


! Print the results to the defined output file
call output()                                         ! tqtec_io.f90


if (verbosity.ge.1) then
    write(*,*) 'tqtec: finished'
    write(*,*) 'Results can be found in ',trim(output_file)
    if (geotherm_file.ne.'') then
        write(*,*) 'Geotherms every ',t_geotherm_output,' Ma can be found in ',trim(geotherm_file)
    endif
endif

end





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

use tqtec, only: verbosity, &
                 nnodes, &
                 dz, &
                 conductivity, &
                 hp, &
                 hf, &
                 nlayers, &
                 layer, &
                 hf_surf, &
                 hf_base, &
                 hp_surf, &
                 hp_dep, &
                 temp_surf, &
                 cond_base, &
                 nhorizons, &
                 temp, &
                 depth, &
                 depth_node

implicit none

! Local variables
integer :: i, j
integer :: ntop, nbot
double precision :: hfhp



if (verbosity.ge.2) then
    write(*,*) 'initialize_thermal_parameters: starting'
endif


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



