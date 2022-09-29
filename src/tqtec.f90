!----
! TQTec (Temperature, Heat Flow, Tectonics)
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
character(len=512) :: temp_file                   ! name of temperature file
character(len=512) :: timing_file                 ! name of tectonic action timing file
integer :: verbosity                              ! name of temperature file

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
! C     W(3) = TEMPERATURE AT BOTTOM NODE + CHANGE IN TEMPERATURE
! C        WITHOUT HEAT PRODUCTION = TEMP AT NODE N+1

! Tectonic events
integer :: nburial                                ! number of burial events                         NBP
double precision, allocatable :: burial_dat(:,:)  ! burial_dat(:,1): start (Ma)                     AN(1)
                                                  ! burial_dat(:,2): duration (Ma)                  AN(2)
                                                  ! burial_dat(:,3): thickness (km)                 AN(3)
                                                  ! burial_dat(:,4): conductivity (W/(m*K))         AN(4)
integer :: nuplift                                ! number of uplift/erosion events                 NUEP
double precision, allocatable :: uplift_dat(:,:)  ! uplift_dat(:,1): start (Ma)                     AN(1)
                                                  ! uplift_dat(:,2): duration (Ma)                  AN(2)
                                                  ! uplift_dat(:,3): thickness (km)                 AN(3)
integer :: nthrust                                ! number of thrusting events                      NTP
double precision, allocatable :: thrust_dat(:,:)  ! thrust_dat(:,1): start (Ma)                     AN(1)
                                                  ! thrust_dat(:,2): upper (1) or lower (2) plate   AN(2)
                                                  ! thrust_dat(:,3): initial base (km)              AZ(1)
                                                  ! thrust_dat(:,4): initial depth (km)             AZ(2)
                                                  ! thrust_dat(:,5): initial thickness (km)         AZ(3)
integer :: nhfvars                                ! number of surface heat flow variations          QSTEP
double precision, allocatable :: hfvar(:,:)       ! (1) start (2) new heat flow                     QVTIME
double precision, allocatable :: bas_grad(:)      ! temperature gradient at the base of the model   BASGRAD
integer, allocatable :: action(:)                 ! burial (1), erosion (2), or thrust (>=3)        P
double precision, allocatable :: bcond(:)         ! boundary condition magnitude                    BCOND

! Results array
double precision, allocatable :: results(:,:,:)   ! temperature and depth for each timstep/depth    r1

end module


!==================================================================================================!
!==================================================================================================!
!==================================================================================================!
!==================================================================================================!
!==================================================================================================!


program main
!----
! Solve for the 1-D transient thermal field defined by boundary conditions and tectonic actions
!----

use tqtec

implicit none

integer :: i, np, ierr
double precision :: cond_surf


! Initialize default model parameters
call initialize_defaults()


! Parse command line
! Matt's note: this is a totally new subroutine for tqtec (which I use in all my other programs)
! that allows better control over user input/output. Here, most of the model I/O is done via a
! control file, so gcmdln() is much simpler, only allowing specification of basic program I/O.
call gcmdln()


if (verbosity.ge.1) then
    write(*,*) 'tqtec: starting'
endif


! Read control file or user input from standard input (formerly INPUT)
call read_model_parameters()


! Calculate model parameters and allocate arrays
if (dz*dble(nnodes).lt.maxval(depth)) then
    write(0,*) 'tqtec: model extent is shallower than deepest horizon'
    call error_exit(1)
endif
nt_total = int(t_total/dt)
nt_geotherm_output = int(t_geotherm_output/dt)
temp_factor = diffusivity*dt/cond_base
r1 = diffusivity*dt/(dz*dz)
dtemp_wo_hp = (hf_surf-hp_surf*hp_dep)*dz/cond_base
allocate(depth_node(nhorizons))
allocate(conductivity(nnodes))
allocate(temp(nnodes))
allocate(hp(nnodes))
allocate(hf(nt_total))
allocate(results(nt_total,2,nhorizons))


! Set up tectonic action timing arrays (formerly: HIST)
call setup_action_arrays()


! Initialize the temperature, heat flow, and heat production at each node (formerly: INIT)
call initialize_thermal_parameters()


! Print model parameters to standard output
if (verbosity.ge.1) then
    call print_model_parameters()
endif


! Print the geotherm to a file (-geotherm flag)
if (temp_file.ne.'') then
    open(unit=12,file=temp_file,status='unknown')
    ! Header contains time step, time since start in Ma, and time until end in Ma
    write(12,'(A,I10,2F10.3)') '> #',0,0.0d0,0.0d0-t_total
    do i = 1,nnodes
        write(12,*) temp(i),dble(i)*dz
    enddo
endif


! Step through time and run the finite difference procedure to calculate the temperature at each
! node and at each time step
istep = 0
do while (istep.lt.nt_total)

    ! Update the adjusted temperature at the base of the model
    temp_base_adj = temp(nnodes) + bas_grad(istep+1)

    ! Calculate the updated temperatures at each node (the main finite difference procedure)
    ! (formerly: MAT and TRID)
    call update_temps(nnodes,ierr)
    if (ierr.ne.0) then
        write(0,*) 'tqtec: error in update_temps() TRID algorithm at step',istep
        stop 1
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

    ! Tectonic action!
    if (action(istep).eq.1) then
        ! Add material to the top of the model (burial)
        call bury() ! (Formerly: BURIAL)

    elseif (action(istep).eq.2) then
        ! Remove material from the top of the model (erosion)
        call erode() ! (Formerly: EROS)

    elseif (action(istep).ge.3) then
        ! Duplicate material and insert it (thrusting)
        if (int(thrust_dat(action(istep)-2,2)).eq.1) then
            ! Track horizons in the hanging wall (upper plate) of the thrust sheet
            call thrust_upperplate() ! (Formerly: THSTUP)
        elseif (int(thrust_dat(action(istep)-2,2)).eq.2) then
            ! Track horizons in the footwall (lower plate) of the thrust sheet
            call thrust_lowerplate() ! (Formerly: THSTLP)
        endif

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
    if (temp_file.ne.'') then
        if (mod(istep,nt_geotherm_output).eq.0) then
            ! Header contains time step, time since start in Ma, and time until end in Ma
            write(12,'(A,I10,2F10.3)') '> #',istep,dble(istep)*dt,dble(istep)*dt-t_total
            do i = 1,nnodes
                write(12,*) temp(i),dble(i)*dz
            enddo
        endif
    endif

enddo


! Close the geotherm file if needed
if (temp_file.ne.'') then
    close(12)
endif


! Print the results to the defined output file
call output()


if (verbosity.ge.1) then
    write(*,*) 'tqtec: finished'
    write(*,*) 'Results can be found in ',trim(output_file)
endif

end







!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!---------------------------------- INPUT SUBROUTINES ---------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!


subroutine initialize_defaults()
!----
! Initialize default model parameters
!----

use tqtec, only: nnodes, &
                 dz, &
                 dt, &
                 nhorizons, &
                 nlayers, &
                 diffusivity, &
                 nburial, &
                 nuplift, &
                 nthrust, &
                 verbosity
implicit none

! Variable = value       ! Value in old tqtec
nnodes = 5000            ! N=1200                 ! number of finite difference spatial nodes
dz = 0.01d0              ! H1=0.05                ! node spacing (km)
dt = 0.001d0             ! K1=0.005               ! time step size (Ma)
nhorizons = 10           ! Hard-coded to 10       ! number of depth horizons to track
nlayers = 0              ! INL                    ! number of layers with different conductivity
diffusivity = 32.0d0     ! D1=32.0                ! thermal diffusivity
nburial = 0              ! NBP                    ! number of burial events
nuplift = 0              ! NUEP                   ! number of uplift/erosion events
nthrust = 0              ! NTP                    ! number of thrust events
verbosity = 0                                     ! program verbosity

return
end subroutine


!--------------------------------------------------------------------------------------------------!


subroutine read_model_parameters()
!----
! Determine how to read input model parameters and run the corresponding input routine
! Depends on value of variable "input_mode":
!     - input_mode="user": interactive input entry
!     - input_mode="file": read input file
!
! Determine how to handle output
!----

use tqtec, only: input_mode, &
                 input_file, &
                 output_file, &
                 verbosity

implicit none

! Local variables
character(len=32) :: reply


if (verbosity.ge.2) then
    write(*,*) 'read_model_parameters: starting'
endif


! Interactive mode (like original version of tqtec)
if (input_mode.eq.'user') then

    ! Ask if user wants to create a new data file
    write(*,*) 'Do you want to manually create a new data file? (Y/N)'
    read(*,*) reply

    if (reply.eq.'y'.or.reply.eq.'Y'.or.reply.eq.'1') then
        write(*,*) 'Name of input file to create?'
        read(*,*) input_file
        call read_interactive()                                 ! Create data file with interactive input
    elseif (reply.eq.'n'.or.reply.eq.'N'.or.reply.eq.'0') then
        write(*,*) 'Name of existing input file?'
        read(*,*) input_file
        call read_input_file()                                  ! Read existing data file
    else
        write(0,*) 'tqtec: could not understand response "',trim(reply),'"'
        stop
    endif


! Read directly from input file
elseif (input_mode.eq.'file') then

    call read_input_file()

else
    write(0,*) 'tqtec: no input mode named "',trim(input_mode),'"'
    stop
endif


! Define an output file if necessary
if (output_file.eq.'') then
    write(*,*) 'Name of output file?'
    read(*,*) output_file
    write(*,*) 'tqtec: creating output file "',trim(output_file),'"'
    write(*,*) 'To create this file automatically, run tqtec -o ',trim(output_file)
endif


if (verbosity.ge.2) then
    write(*,*) 'read_model_parameters: finished'
endif

return
end subroutine


!--------------------------------------------------------------------------------------------------!


subroutine read_interactive()
!----
! Manually enter model parameters and tectonic events
!----

use tqtec, only: input_file, &
                 t_total, &
                 t_geotherm_output, &
                 temp_surf, &
                 hf_surf, &
                 hp_surf, &
                 hp_dep, &
                 cond_base, &
                 nlayers, &
                 layer, &
                 nhorizons, &
                 depth, &
                 nburial, &
                 burial_dat, &
                 nuplift, &
                 uplift_dat, &
                 nthrust, &
                 thrust_dat, &
                 nhfvars, &
                 hfvar

implicit none

! Local variables
integer :: i, j
character(len=32) :: reply, fmt_string


! Open the input file so model parameters can be saved to it
open(unit=9,file=input_file,status='unknown')


! Write to fixed format input file
write(9,*) trim(input_file)


! Model timing
write(*,*) 'Total time for model? (Ma)'
read(*,*) t_total
write(*,*) 'Time interval between geotherm outputs (if -geotherm flag is used)? (Ma)'
read(*,*) t_geotherm_output

! Model boundary conditions
write(*,*) 'Temperature at upper surface boundary? (C)'
read(*,*) temp_surf
write(*,*) 'Surface heat flow? (mW/m^2)'
read(*,*) hf_surf
write(*,*) 'Initial (basement) thermal conductivity? (W/(m*K))'
read(*,*) cond_base
write(*,*) 'Surface heat production? (uW/m^3)'
read(*,*) hp_surf
write(*,*) 'Heat production depth? (km)'
read(*,*) hp_dep

! Write to fixed format input file
write(9,1001) t_total, t_geotherm_output, temp_surf, hf_surf, cond_base, hp_surf, hp_dep
1001 format(2F10.0,5F10.4)


! Variations in thermal conductivity
write(*,*) 'Do you want to account for variations in thermal conductivity at the start ',&
           'of the model? (y/n)'
read(*,*) reply
if (reply.eq.'y'.or.reply.eq.'Y'.or.reply.eq.'1') then
    write(*,*) 'Number of layers to input conductivity for?'
    read(*,*) nlayers
    do i = 1,nlayers
        write(*,*) 'Depth of top of layer',i,'? (km)'
        read(*,*) layer(i,1)
        write(*,*) 'Thickness of layer',i,'? (km)'
        read(*,*) layer(i,2)
        write(*,*) 'Conductivity of layer',i,'? (W/(m*K))'
        read(*,*) layer(i,3)
    enddo
endif

! Write to fixed format input file
write(9,'(I10)') nlayers
if (nlayers.gt.0) then
    do i = 1,nlayers
        write(9,1003) (layer(i,j),j=1,3)
    enddo
endif
1003 format(3F10.4)


! Tracked horizon depths
write(*,*) 'Number of horizons to track? (Press <return> to use default: 10)'
read(*,'(A)') reply
if (reply.ne.'') then
    read(reply,*) nhorizons
endif
if (allocated(depth)) then
    deallocate(depth)
endif
allocate(depth(nhorizons))
do i = 1,nhorizons
    write(*,*) 'Initial depth of point/horizon',i,' (of ',nhorizons,')? (km)'
    read(*,*) depth(i)
enddo

! Write to fixed format input file
write(fmt_string,'("(",I5,"F8.4",")")') nhorizons
write(9,fmt_string) (depth(i),i=1,nhorizons)


! Burial events
write(*,*) 'Number of burial periods?'
read(*,*) nburial
if (allocated(burial_dat)) then
    deallocate(burial_dat)
endif
allocate(burial_dat(nburial,4))
do i = 1,nburial
    write(*,*) 'Beginning of burial period',i,'? (Ma after start)'
    read(*,*) burial_dat(i,1)
    write(*,*) 'Duration of burial period',i,'? (Ma)'
    read(*,*) burial_dat(i,2)
    write(*,*) 'Total burial during episode',i,'? (km)'
    read(*,*) burial_dat(i,3)
    write(*,*) 'Thermal conductivity of sediments in burial episode',i,'? (W/(m*K))'
    read(*,*) burial_dat(i,4)
enddo

! Write to fixed format input file
write(9,'(I10)') nburial
if (nburial.gt.0) then
    do i = 1,nburial
        write(9,'(4F10.4)') (burial_dat(i,j),j=1,4)
    enddo
endif


! Uplift events
write(*,*) 'Number of uplift/erosion periods?'
read(*,*) nuplift
if (allocated(uplift_dat)) then
    deallocate(uplift_dat)
endif
allocate(uplift_dat(nuplift,3))
do i = 1,nuplift
    write(*,*) 'Beginning of uplift period',i,'? (Ma after start)'
    read(*,*) uplift_dat(i,1)
    write(*,*) 'Duration of uplift period',i,'? (Ma)'
    read(*,*) uplift_dat(i,2)
    write(*,*) 'Total uplift during episode',i,'? (km)'
    read(*,*) uplift_dat(i,3)
enddo

! Write to fixed format input file
write(9,'(I10)') nuplift
if (nuplift.gt.0) then
    do i = 1,nuplift
        write(9,'(3F10.4)') (uplift_dat(i,j),j=1,3)
    enddo
endif


! Thrust events
write(*,*) 'Number of thrust periods?'
read(*,*) nthrust
if (allocated(thrust_dat)) then
    deallocate(thrust_dat)
endif
allocate(thrust_dat(nthrust,5))
do i = 1,nthrust
    write(*,*) 'Time of thrust period',i,'? (Ma after start)'
    read(*,*) thrust_dat(i,1)
    write(*,*) 'Points in upper(1) or lower(2) plate during thrust period',i,'?'
    read(*,*) thrust_dat(i,2)
    write(*,*) 'Initial base of thrust during episode',i,'? (km)'
    read(*,*) thrust_dat(i,3)
    write(*,*) 'Initial depth of thrust during episode',i,'? (km)'
    read(*,*) thrust_dat(i,4)
    write(*,*) 'Initial thickness of thrust during episode',i,'? (km)'
    read(*,*) thrust_dat(i,5)
enddo

! Write to fixed format input file
write(9,'(I10)') nthrust
if (nthrust.gt.0) then
    do i = 1,nthrust
        write(9,'(5F10.4)') (thrust_dat(i,j),j=1,5)
    enddo
endif


! Variations in heat flow
write(*,*) 'Number of heat flow variations?'
read(*,*) nhfvars
if (allocated(hfvar)) then
    deallocate(hfvar)
endif
allocate(hfvar(nhfvars,2))
do i = 1,nhfvars
    write(*,*) 'Time of heat flow value change',i,'? (Ma after start)'
    read(*,*) hfvar(i,1)
    write(*,*) 'Value of heat flow at change',i,'?'
    read(*,*) hfvar(i,2)
enddo

! Write to fixed format input file
write(9,'(I10)') nhfvars
if (nhfvars.gt.0) then
    do i = 1,nhfvars
        write(9,'(2F10.4)') (hfvar(i,j),j=1,2)
    enddo
endif


write(*,*) 'tqtec: input file "',trim(input_file),'" has been created'
write(*,*) 'To re-use this file, run tqtec -f ',trim(input_file)


return
end subroutine


!--------------------------------------------------------------------------------------------------!


subroutine read_input_file()
!----
! Determine whether to read fixed format (original) or free format (new) input file
!----

use tqtec, only: input_file

implicit none

! Local variables
integer :: ios
character(len=512) :: input_line
logical :: ex



! Check to make sure input file exists
inquire(file=input_file,exist=ex)
if (.not.ex) then
    write(0,*) 'tqtec: could not find input file "',trim(input_file),'"'
    stop
endif

! Check the input file format
! Old version is fixed format
! New version has VAR=VALUE
open(unit=8,file=input_file,iostat=ios)
if (ios.ne.0) then
    write(0,*) 'tqtec: something went wrong trying to open input file "',trim(input_file),'"'
    stop
endif
read(8,'(A)') input_line
close(8)
if (index(input_line,'=').ne.0.or.input_line(1:1).eq.'#') then
    call read_input_file_new()
else
    call read_input_file_old()
endif


return
end subroutine


!--------------------------------------------------------------------------------------------------!


subroutine read_input_file_old()
!----
! Read tqtec input file in original fixed format, e.g.:
!
! tqtec.in
!         50         5    0.0000   30.0000    3.0000    0.0000
!          0
!   2.0000  4.0000  6.0000  8.0000 10.0000 12.0000 14.0000 16.0000 18.0000 20.0000
!          1
!    10.0000   10.0000    5.0000    2.0000
!          1
!    20.0000   20.0000   10.0000
!          1
!    40.0000         1   25.0000    0.0000  25.0000
!----

use tqtec, only: input_file, &
                 verbosity, &
                 t_total, &
                 t_geotherm_output, &
                 temp_surf, &
                 hf_surf, &
                 hp_surf, &
                 hp_dep, &
                 cond_base, &
                 nlayers, &
                 layer, &
                 nhorizons, &
                 depth, &
                 nburial, &
                 burial_dat, &
                 nuplift, &
                 uplift_dat, &
                 nthrust, &
                 thrust_dat, &
                 nhfvars, &
                 hfvar

implicit none

! Local variables
integer :: i, j, ios
character(len=512) :: input_line
logical :: inWhitespace


if (verbosity.ge.3) then
    write(*,*) '    read_input_file_old: starting'
endif


! Open the input file for reading in old fixed format
open(unit=8,file=input_file,iostat=ios)
if (ios.ne.0) then
    write(0,*) 'tqtec: something went wrong trying to open input file "',trim(input_file),'"'
    stop
endif


ios = 0

! First line contains file name
read(8,'(A)') input_line


! Second line contains model parameters
! READ (8,110) Q1,M1,W(1),G1,C1,A1,B1
read(8,'(A)') input_line
read(input_line,*,iostat=ios) t_total, t_geotherm_output, temp_surf, hf_surf, cond_base, hp_surf, &
                   hp_dep
if (ios.ne.0) then
    read(input_line,*,iostat=ios) t_total, t_geotherm_output, temp_surf, hf_surf, cond_base, hp_surf
endif
if (ios.ne.0) then
    write(0,*) 'tqtec: something went wrong trying to read model parameters'
    write(0,*) 'Offending line:'
    write(0,*) trim(input_line)
    stop 1
endif


! Read material layers
! READ (8,150) INL
read(8,*) nlayers
if (allocated(layer)) then
    deallocate(layer)
endif
allocate(layer(nlayers,3))
do i = 1,nlayers
    read(8,*) (layer(i,j),j=1,3)
enddo


! Read horizon depths
! Any number of horizons can be listed here, so reset nhorizons and deallocate depth array
nhorizons = 0
if (allocated(depth)) then
    deallocate(depth)
endif
read(8,'(A)') input_line
! Parse the input line for the number of depth horizons
i = 1
inWhitespace = .true.
do while (i.le.len_trim(input_line))
    if (input_line(i:i).eq.' ') then
        inWhitespace = .true.
    else
        if (inWhitespace) then
            nhorizons = nhorizons + 1
        endif
        inWhitespace = .false.
    endif
    i = i + 1
enddo
! Reallocate depth array and read depths
allocate(depth(nhorizons))
read(input_line,*) (depth(i),i=1,nhorizons)


! Read burial episodes
read(8,*,end=1001,iostat=ios) nburial
if (allocated(burial_dat)) then
    deallocate(burial_dat)
endif
allocate(burial_dat(nburial,4))
do i = 1,nburial
    read(8,'(A)',end=1101,iostat=ios) input_line
    read(input_line,*,end=1201,iostat=ios) (burial_dat(i,j),j=1,4)
enddo


! Read uplift/erosion episodes
read(8,*,end=1002,iostat=ios) nuplift
if (allocated(uplift_dat)) then
    deallocate(uplift_dat)
endif
allocate(uplift_dat(nuplift,4))
do i = 1,nuplift
    read(8,'(A)',end=1102,iostat=ios) input_line
    read(input_line,*,end=1202,iostat=ios) (uplift_dat(i,j),j=1,3)
enddo


! Read thrust episodes
read(8,*,end=1003,iostat=ios) nthrust
if (allocated(thrust_dat)) then
    deallocate(thrust_dat)
endif
allocate(thrust_dat(nthrust,5))
do i = 1,nthrust
    read(8,'(A)',end=1103,iostat=ios) input_line
    read(input_line,*,end=1203,iostat=ios) (thrust_dat(i,j),j=1,5)
enddo


! Read basal heat flow variations
read(8,*,end=1004,iostat=ios) nhfvars
if (allocated(hfvar)) then
    deallocate(hfvar)
endif
allocate(hfvar(nhfvars,2))
do i = 1,nhfvars
    read(8,'(A)',end=1104,iostat=ios) input_line
    read(input_line,*,end=1204,iostat=ios) (hfvar(i,j),j=1,2)
enddo


! Warning messages if unable to find tectonic events
1001 if (ios.ne.0) then
    write(0,*) 'tqtec: could not find any burial event settings in input file'
endif
1002 if (ios.ne.0) then
    write(0,*) 'tqtec: could not find any uplift event settings in input file'
endif
1003 if (ios.ne.0) then
    write(0,*) 'tqtec: could not find any thrust event settings in input file'
endif
1004 if (ios.ne.0) then
    write(0,*) 'tqtec: could not find any heat flow variation settings in input file'
endif
ios = 0

! Errors if unable to read number of specified tectonic events
1101 if (ios.ne.0) then
    write(0,*) 'tqtec: input file only specified',i-1,' of',nburial,' burial events'
    stop 1
endif
1102 if (ios.ne.0) then
    write(0,*) 'tqtec: input file only specified',i-1,' of',nuplift,' uplift events'
    stop 1
endif
1103 if (ios.ne.0) then
    write(0,*) 'tqtec: input file only specified',i-1,' of',nthrust,' thrust events'
    stop 1
endif
1104 if (ios.ne.0) then
    write(0,*) 'tqtec: input file only specified',i-1,' of',nhfvars,' heat flow variations'
    stop 1
endif

! Errors if tectonic events lines are too short
1201 if (ios.ne.0) then
    write(0,*) 'tqtec: could not parse:'
    write(0,*) '"',trim(input_line),'"'
    write(0,*) 'As: TSTART  TDURATION  BURIAL  CONDUCTIVITY'
    stop 1
endif
1202 if (ios.ne.0) then
    write(0,*) 'tqtec: could not parse:'
    write(0,*) '"',trim(input_line),'"'
    write(0,*) 'As: TSTART  TDURATION  UPLIFT'
    stop 1
endif
1203 if (ios.ne.0) then
    write(0,*) 'tqtec: could not parse:'
    write(0,*) '"',trim(input_line),'"'
    write(0,*) 'As: TSTART  1/0  BASE  DEPTH  THICKNESS'
    stop 1
endif
1204 if (ios.ne.0) then
    write(0,*) 'tqtec: could not parse:'
    write(0,*) '"',trim(input_line),'"'
    write(0,*) 'As: TSTART  HEAT_FLOW'
    stop 1
endif


close(8)

if (verbosity.ge.3) then
    write(*,*) '    read_input_file_old: finished'
endif

return
end subroutine


!--------------------------------------------------------------------------------------------------!


subroutine read_input_file_new()
!----
! Read tqtec input file in more flexible format, e.g.:
!
! T_TOTAL=50
! T_GEOTHERM_OUTPUT=5
! TEMP_SURF=0
! HF_SURF=30
! COND_BASE=3
! HP_SURF=0
! HP_DEP=0
! NLAYERS=0
! NHORIZONS=10
! 2 4 6 8 10 12 14 16 18 20
! NBURIAL=1
! 10 10 5 2
! NUPLIFT=2
! 20 20 10
! 45 2 1
! NTHRUST=1
! 40 1 25 0 25
! NHFVARS=1
! 15 10
! NNODES=
! DZ=
! MAX_DEPTH=
!----

use tqtec

implicit none

! Local variables
integer :: i, j, ios, iend
character(len=32) :: var, value
character(len=512) :: input_line
double precision :: max_depth
logical :: isMaxDepthDefined



ios = 0
iend = 0
isMaxDepthDefined = .false.


! Open the input file for reading in free format
open(unit=8,file=input_file,iostat=ios)
if (ios.ne.0) then
    write(0,*) 'tqtec: something went wrong trying to open input file "',trim(input_file),'"'
    stop
endif


! Initialize required variables so they can later be checked
t_total = -1.0d99
t_geotherm_output = -1.0d99
temp_surf = -1.0d99
hf_surf = -1.0d99
cond_base = -1.0d99
hp_surf = 0.0d0
hp_dep = 0.0d0
nhorizons = 0


! Read the file in flexible format
do while (iend.eq.0)
    read(8,'(A)',end=3451,iostat=iend) input_line

    ! Lines that start with # are comments
    if (input_line(1:1).eq.'#') then
        cycle
    endif

    ! All variables are in the format VAR=VALUE
    i = index(input_line,'=')
    input_line(i:i) = ' '
    read(input_line,*,iostat=ios) var, value
    if (ios.ne.0) then
        write(0,*) 'tqtec: something went wrong trying to read "',trim(input_line),'"'
        stop
    endif

    ! Big if statement to handle all cases of VAR definitions
    if (var.eq.'T_TOTAL'.or.var.eq.'t_total') then
        read(value,*) t_total
    elseif (var.eq.'T_GEOTHERM_OUTPUT'.or.var.eq.'t_geotherm_output') then
        read(value,*) t_geotherm_output
    elseif (var.eq.'TEMP_SURF'.or.var.eq.'temp_surf') then
        read(value,*) temp_surf
    elseif (var.eq.'HF_SURF'.or.var.eq.'hf_surf') then
        read(value,*) hf_surf
    elseif (var.eq.'COND_BASE'.or.var.eq.'cond_base') then
        read(value,*) cond_base
    elseif (var.eq.'HP_SURF'.or.var.eq.'hp_surf') then
        read(value,*) hp_surf
    elseif (var.eq.'HP_DEP'.or.var.eq.'hp_dep') then
        read(value,*) hp_dep
    elseif (var.eq.'NLAYERS'.or.var.eq.'nlayers') then
        read(value,*) nlayers
        if (nlayers.gt.0) then
            if (allocated(layer)) then
                deallocate(layer)
            endif
            allocate(layer(nlayers,3))
            do i = 1,nlayers
                read(8,*) (layer(i,j),j=1,3)                ! top thickness conductivity
            enddo
        endif
    elseif (var.eq.'NHORIZONS'.or.var.eq.'nhorizons') then
        read(value,*) nhorizons
        if (nhorizons.gt.0) then
            if (allocated(depth)) then
                deallocate(depth)
            endif
            allocate(depth(nhorizons))
            read(8,*) (depth(i),i=1,nhorizons)              ! depth
        endif
    elseif (var.eq.'NBURIAL'.or.var.eq.'nburial') then
        read(value,*) nburial
        if (nburial.gt.0) then
            if (allocated(burial_dat)) then
                deallocate(burial_dat)
            endif
            allocate(burial_dat(nburial,4))
            do i = 1,nburial
                read(8,*) (burial_dat(i,j),j=1,4)           ! start duration thickness conductivity
            enddo
        endif
    elseif (var.eq.'NUPLIFT'.or.var.eq.'nuplift') then
        read(value,*) nuplift
        if (nuplift.gt.0) then
            if (allocated(uplift_dat)) then
                deallocate(uplift_dat)
            endif
            allocate(uplift_dat(nuplift,3))
            do i = 1,nuplift
                read(8,*) (uplift_dat(i,j),j=1,3)           ! start duration thickness
            enddo
        endif
    elseif (var.eq.'NTHRUST'.or.var.eq.'nthrust') then
        read(value,*) nthrust
        if (nthrust.gt.0) then
            if (allocated(thrust_dat)) then
                deallocate(thrust_dat)
            endif
            allocate(thrust_dat(nthrust,5))
            do i = 1,nthrust
                read(8,*) (thrust_dat(i,j),j=1,5)           ! start upper/lower thick_init depth thick_final
            enddo
        endif
    elseif (var.eq.'NHFVARS'.or.var.eq.'nhfvars') then
        read(value,*) nhfvars
        if (nhfvars.gt.0) then
            if (allocated(hfvar)) then
                deallocate(hfvar)
            endif
            allocate(hfvar(nhfvars,2))
            do i = 1,nhfvars
                read(8,*) (hfvar(i,j),j=1,2)                ! start heat_flow
            enddo
        endif
    elseif (var.eq.'NNODES'.or.var.eq.'nnodes') then
        read(value,*) nnodes
    elseif (var.eq.'DZ'.or.var.eq.'dz') then
        read(value,*) dz
        if (isMaxDepthDefined) then
            nnodes = int(max_depth/dz)
        endif
    elseif (var.eq.'MAX_DEPTH'.or.var.eq.'max_depth'.or.var.eq.'DEPTH_MAX'.or.var.eq.'depth_max') then
        read(value,*) max_depth
        nnodes = int(max_depth/dz)
        isMaxDepthDefined = .true.
    else
        write(0,*) 'tqtec: no variable option named "',trim(var),'"'
        stop
    endif

    ! Reached the end of the file, exit
    3451 if (iend.ne.0) then
        exit
    endif
enddo


! Check that necessary variables have been defined
if (t_total.lt.0.0) then
    write(0,*) 'tqtec: t_total has not been defined'
    stop
endif
if (temp_surf.lt.0.0) then
    write(0,*) 'tqtec: temp_surf has not been defined'
    stop
endif
if (hf_surf.lt.0.0) then
    write(0,*) 'tqtec: hf_surf has not been defined'
    stop
endif
if (cond_base.lt.0.0) then
    write(0,*) 'tqtec: cond_base has not been defined'
    stop
endif
if (nhorizons.le.0) then
    write(0,*) 'tqtec: no horizons have been defined'
    stop
endif
if (temp_surf.lt.0.0) then
    write(0,*) 'tqtec: temp_surf has not been defined'
    stop
endif
if (nburial.eq.0.and.nuplift.eq.0.and.nthrust.eq.0.and.nhfvars.eq.0) then
    write(0,*) 'tqtec: no tectonic actions have been defined'
    stop
endif

close(8)

return
end subroutine






!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!---------------------------------- MODEL PREPARATION ---------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!


subroutine setup_action_arrays()
!----
! Define arrays to control burial, erosion, and thrusting events
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
                 nhfvars, &
                 hfvar, &
                 bas_grad, &
                 action, &
                 bcond

implicit none


! Local variables
integer :: i, ct
integer :: j, jbeg, jend
integer :: nstart, nduration, nthick
double precision :: rate, arg
double precision :: hf_surf_var(nt_total)
integer :: intqt(nhfvars)


if (verbosity.ge.2) then
    write(*,*) 'setup_action_arrays: starting'
endif


! Allocate memory to tectonic action arrays
if (allocated(action)) then
    deallocate(action)
endif
if (allocated(bcond)) then
    deallocate(bcond)
endif
allocate(action(nt_total))
allocate(bcond(nt_total))
allocate(bas_grad(nt_total))


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
            action(j) = 1
            bcond(j) = burial_dat(i,4)
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
            action(j) = 2
            ct = ct + 1
        endif
    enddo
enddo


! Thrust events
do i = 1,nthrust
    nstart = int(thrust_dat(i,1)/dt)              ! Starting timestep
    action(nstart) = i+2
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


subroutine print_model_parameters()
!----
! Print salient model parameters to standard output (useful for debugging)
!----

use tqtec

implicit none

! Local variables
integer :: i, j

write(*,*) 'Nodes'
write(*,2002) 'nnodes:           ',nnodes
write(*,2001) 'dz:               ',dz,'km'
write(*,2001) 'max_depth:        ',dble(nnodes)*dz,'km'
write(*,*) 'Timing'
write(*,2001) 't_total:          ',t_total,'Ma'
write(*,2002) 'nt_total:         ',nt_total
write(*,2001) 't_geotherm_output:',t_geotherm_output,'Ma'
write(*,2001) 'dt:               ',dt,'Ma'
write(*,*) 'Boundary conditions'
write(*,2001) 'temp_surf:        ',temp_surf,'C'
write(*,2001) 'hf_surf:          ',hf_surf,'mW/m^2'
write(*,2001) 'hp_surf:          ',hp_surf,'uW/m^3'
write(*,2001) 'hp_dep:           ',hp_dep,'km'
write(*,*) 'Material properties'
write(*,2001) 'cond_base:        ',cond_base,'W/(m*K)'
write(*,2001) 'diffusivity:      ',diffusivity,'km^2/Ma'
write(*,2002) 'nlayers:          ',nlayers
if (nlayers.gt.0) then
    write(*,'(5X,3A14)') 'top(km)', 'thick(km)', 'cond(W/(m*K))'
    do i = 1,nlayers
        write(*,'(5X,3F14.3)') layer(i,1),layer(i,2),layer(i,3)
    enddo
endif
write(*,*) 'Tracked horizons'
write(*,2002) 'nhorizons:        ',nhorizons
write(*,'(5X,4A14)') 'depth(km)','depth_node'
do i = 1,nhorizons
    write(*,'(5X,F14.3,I14)') depth(i),depth_node(i)
enddo
write(*,*) 'Tectonic actions'
write(*,2002) 'nburial:          ',nburial
if (nburial.gt.0) then
    write(*,'(5X,4A14)') 'start(Ma)', 'duration(Ma)', 'thickness(km)', 'cond(W/(m*K))'
    do i = 1,nburial
        write(*,'(5X,4F14.3)') (burial_dat(i,j),j=1,4)
    enddo
endif
write(*,2002) 'nuplift:          ',nuplift
if (nuplift.gt.0) then
    write(*,'(5X,3A14)') 'start(Ma)', 'duration(Ma)', 'thickness(km)'
    do i = 1,nuplift
        write(*,'(5X,3F14.3)') (uplift_dat(i,j),j=1,3)
    enddo
endif
write(*,2002) 'nthrust:          ',nthrust
if (nthrust.gt.0) then
    write(*,'(5X,2A14,2X,3A14)') 'start(Ma)', 'upper/lower', 'thick_init(km)', 'dep_base(km)', &
                                 'thick_end(km)'
    do i = 1,nthrust
        write(*,'(5X,2F14.3,2X,3F14.3)') (thrust_dat(i,j),j=1,5)
    enddo
endif

2001 format(5X,A18,F10.3,X,A)
2002 format(5X,A18,I10,X,A)

! write(0,*) 'istep:         ',istep
! write(0,*) 'r1:            ',r1
! write(0,*) 'nt_geotherm_output:     ',nt_geotherm_output
! write(0,*) 'hf_base:       ',hf_base
! write(0,*) 'dtemp_wo_hp:   ',dtemp_wo_hp
! write(0,*) 'dtemp_wo_hp:   ',dtemp_wo_hp
! write(0,*) 'temp_factor:   ',temp_factor
! write(0,*) 'temp_base_adj:  ',temp_base_adj
! write(0,*) 'nhfvars:        ',nhfvars
! write(0,*) 'hfvar(:,1):     ',hfvar(:,1)
! write(0,*) 'hfvar(:,2):     ',hfvar(:,2)

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



!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!----------------------------------- TECTONIC ACTIONS ---------------------------------------------!
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
                 bcond

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
conductivity(1) = bcond(istep)

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
                 thrust_dat, &
                 action

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
k = action(istep) - 2

! Save thrust sheet parameters
thick_init = int(thrust_dat(k,3)/dz)  ! Thickness prior to thrusting, in nodes
thrust_dep = int(thrust_dat(k,4)/dz)  ! Depth of emplacement, in nodes
thick_end = int(thrust_dat(k,5)/dz)   ! Final thickness of thrust sheet, in nodes
if (thick_init.lt.thick_end) then
    write(0,*) 'thrust_upperplate: final thickness must be less than or equal to initial thickness'
    stop 1
endif

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
                 thrust_dat, &
                 action

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
k = action(istep) - 2

! Save thrust sheet parameters
thick_init = int(thrust_dat(k,3)/dz)  ! Thickness prior to thrusting, in nodes
thrust_dep = int(thrust_dat(k,4)/dz)  ! Depth of emplacement, in nodes
thick_end = int(thrust_dat(k,5)/dz)   ! Final thickness of thrust sheet, in nodes
if (thick_init.lt.thick_end) then
    write(0,*) 'thrust_lowerplate: final thickness must be less than or equal to initial thickness'
endif

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
!--------------------------------------------------------------------------------------------------!
!-------------------------------------- OUTPUTS ---------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!


subroutine output()
!----
! Print model results to a file
!----

use tqtec

implicit none

! Local variables
integer :: i, j, k


! Open the output file
open(unit=7,file=output_file,status='unknown')


! Write the results in the format originally specified by Kevin Furlong

! Output file name
write(7,*) trim(output_file)

! Model parameters
write(7,110) dz
write(7,110) dt
write(7,110) hp_surf
write(7,110) hp_dep
write(7,110) t_total
write(7,110) diffusivity
write(7,110) temp_factor
write(7,'(I10)') nhorizons
write(7,110) 0.0 ! II(9)
write(7,110) 0.0 ! II(10)
110 format(F7.3)

! Heat flow
do j = 2,nt_total,2
    write(7,115) hf(j)
enddo
115 format(F6.2)

! Temperature and depth of tracked horizons
do k = 1,nhorizons
    do j = 1,2
        do i = 2,nt_total,2
            write(7,120) results(i,j,k)
        enddo
    enddo
enddo
120 format(F7.1)

! Horizon depths
do i = 1,nhorizons
    write(7,130) depth_node(i)
enddo
! 130 format(F11.4)
130 format(I11)

close(7)

return
end subroutine


!--------------------------------------------------------------------------------------------------!


subroutine gcmdln()
!----
! Parse tqtec command line arguments defining input/output modes and files
!----

use tqtec, only: input_mode, &
                 input_file, &
                 output_file, &
                 temp_file, &
                 timing_file, &
                 verbosity, &
                 nnodes, &
                 dz, &
                 dt

implicit none

! Local variables
character(len=512) arg
integer :: i, j, ios, narg


! Initialize control variables
ios = 0

! Initialize defaults
input_file = ''
input_mode = 'user'
output_file = ''
temp_file = ''
timing_file = ''


narg = command_argument_count()
if (narg.eq.0) then
    call usage('')
endif


i = 1
do while (i.le.narg)

    call get_command_argument(i,arg)

    if (arg.eq.'-f') then
        input_mode = 'file'
        i = i + 1
        call get_command_argument(i,input_file,status=ios)

    elseif (arg.eq.'-i'.or.arg.eq.'-interactive') then
        input_mode = 'user'

    elseif (arg.eq.'-o') then
        i = i + 1
        call get_command_argument(i,output_file,status=ios)

    elseif (arg.eq.'-geotherm') then
        i = i + 1
        call get_command_argument(i,temp_file,status=ios)

    elseif (arg.eq.'-timing') then
        i = i + 1
        call get_command_argument(i,timing_file,status=ios)

    elseif (arg.eq.'-v'.or.arg.eq.'-verbosity') then
        i = i + 1
        call get_command_argument(i,arg,status=ios)
        read(arg,*) verbosity

    elseif (arg.eq.'-h') then
        call usage('')

    elseif (arg(1:7).eq.'NNODES=') then
        j = index(arg,'=')
        arg(1:j) = ' '
        read(arg,*,iostat=ios) nnodes
    elseif (arg(1:3).eq.'DZ=') then
        j = index(arg,'=')
        arg(1:j) = ' '
        read(arg,*,iostat=ios) dz
    elseif (arg(1:3).eq.'DT=') then
        j = index(arg,'=')
        arg(1:j) = ' '
        read(arg,*,iostat=ios) dt

    else
        call usage('tqtec: no option '//trim(arg))
    endif

    if (ios.ne.0) then
        call usage('tqtec: error parsing "'//trim(arg)//'" flag arguments')
    endif

    i = i + 1
enddo

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine usage(str)
implicit none
character(len=*) :: str
if (str.ne.'') then
    write(0,*) trim(str)
    write(0,*)
endif
write(0,*) 'Usage: tqtec -i|-f INPUT_FILE  [-o OUTPUT_FILE] [-geotherm TEMP_FILE] [-timing TIMING_FILE]'
write(0,*)
write(0,*) '-i[nteractive]        Interactively defined model parameters'
write(0,*) '-f INPUT_FILE         Input model parameter file'
write(0,*) '-o OUTPUT_FILE        Output temperature-depth-time file for specified horizons'
write(0,*) '-geotherm TEMP_FILE   Geotherms (output frequency defined in INPUT_FILE)'
write(0,*) '-timing Timing_FILE   Timing of tectonic actions'
write(0,*) '-v VERBOSITY          Verbosity level'
write(0,*)
call error_exit(1)
stop
end subroutine
