!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!---------------------------------- INPUT SUBROUTINES ---------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!




subroutine gcmdln()
!----
! Parse tqtec command line arguments defining input/output modes and files
!----

use tqtec, only: input_file, &
                 input_mode, &
                 output_file, &
                 geotherm_file, &
                 timing_file, &
                 verbosity, &
                 nnodes, &
                 dz, &
                 dt

implicit none

! Local variables
character(len=512) :: arg
character(len=8) :: exec_name
integer :: i, j, ios, narg


! Initialize control variables
ios = 0


! Initialize default values
input_file = ''
input_mode = 'file'
output_file = ''
geotherm_file = ''
timing_file = ''
#ifdef COMPILE_TQTEC
    exec_name = 'tqtec'
#elif COMPILE_TQCLIM
    exec_name = 'tqclim'
#else
    exec_name = 'NAME_ERR'
#endif


! Count arguments, then exit with usage statement if no arguments
narg = command_argument_count()
if (narg.eq.0) then
    call usage('')
endif



! Parse command line arguments
i = 1
do while (i.le.narg)


    call get_command_argument(i,arg)


    ! Input file
    if (arg.eq.'-f') then
        input_mode = 'file'
        i = i + 1
        call get_command_argument(i,input_file,status=ios)


    ! Print input file format
    elseif (arg.eq.'-f:d'.or.arg.eq.'-f:details'.or.arg.eq.'-f:example') then
        call print_input_file_details()


    ! Input interactively (not recommended)
    elseif (arg.eq.'-i'.or.arg.eq.'-interactive') then
        input_mode = 'user'
        write(0,*) 'WARNING: Interactive mode is no longer supported and may not work properly'
        write(0,*) 'We recommend that you use an input file (-f)'


    ! Output file
    elseif (arg.eq.'-o') then
        i = i + 1
        call get_command_argument(i,output_file,status=ios)


    ! Geotherm file
    elseif (arg.eq.'-geotherm') then
        i = i + 1
        call get_command_argument(i,geotherm_file,status=ios)


    ! Timing action file
    elseif (arg.eq.'-timing') then
        i = i + 1
        call get_command_argument(i,timing_file,status=ios)


    ! Verbosity level
    elseif (arg.eq.'-v'.or.arg.eq.'-verbosity') then
        i = i + 1
        call get_command_argument(i,arg,status=ios)
        read(arg,*) verbosity


    ! Print usage
    elseif (arg.eq.'-h') then
        call usage('')


    ! Set model parameters
    elseif (arg(1:7).eq.'NNODES=') then
        j = index(arg,'=')
        arg(1:j) = ' '
        read(arg,*,iostat=ios) nnodes
        write(0,*) 'WARNING: We recommend setting NNODES in the input file instead of the command line'
    elseif (arg(1:3).eq.'DZ=') then
        j = index(arg,'=')
        arg(1:j) = ' '
        read(arg,*,iostat=ios) dz
        write(0,*) 'WARNING: We recommend setting DZ in the input file instead of the command line'
    elseif (arg(1:3).eq.'DT=') then
        j = index(arg,'=')
        arg(1:j) = ' '
        read(arg,*,iostat=ios) dt
        write(0,*) 'WARNING: We recommend setting DT in the input file instead of the command line'


    ! Could not find this option...
    else
        call usage(trim(exec_name)//': no option '//trim(arg))
    endif


    ! Ran into an error...
    if (ios.ne.0) then
        call usage(trim(exec_name)//': error parsing "'//trim(arg)//'" flag arguments')
    endif


    i = i + 1

enddo

return
end subroutine



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
#ifdef COMPILE_TQTEC
dz = 0.01d0              ! H1=0.05                ! node spacing (km)
dt = 0.001d0             ! K1=0.005               ! time step size (Ma)
#elif COMPILE_TQCLIM
dz = 1.0d0                                        ! node spacing (m)
dt = 1.0d0                                        ! time step size (yr)
#endif
nhorizons = 10           ! Hard-coded to 10       ! number of depth horizons to track
nlayers = 0              ! INL                    ! number of layers with different conductivity
diffusivity = 32.0d0     ! D1=32.0                ! thermal diffusivity (m^2/yr = km^2/Ma)
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
character(len=8) :: exec_name


if (verbosity.ge.2) then
    write(*,*) 'read_model_parameters: starting'
endif


! Set name of executable
#ifdef COMPILE_TQTEC
    exec_name = 'tqtec'
#elif COMPILE_TQCLIM
    exec_name = 'tqclim'
#else
    exec_name = 'NAME_ERR'
#endif




if (input_mode.eq.'user') then  ! Interactive mode (like original tqtec)

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
        write(0,*) trim(exec_name)//': could not understand response "',trim(reply),'"...'
        write(0,*) 'Exiting '//trim(exec_name)
        call error_exit(1)
    endif



elseif (input_mode.eq.'file') then  ! Input file mode (for batch processing)

    call read_input_file()



else

    write(0,*) trim(exec_name)//': no input mode named "',trim(input_mode),'"'
    call error_exit(1)

endif


! Define an output file if necessary
if (output_file.eq.'') then
    write(*,*) 'Name of output file?'
    read(*,*) output_file
    write(*,*) trim(exec_name)//': creating output file "',trim(output_file),'"'
    write(*,*) 'To create this file automatically, run '//trim(exec_name)//' -o ',trim(output_file)
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



write(0,*) 'WARNING: reading model parameters interactively is no longer supported'
write(0,*) 'You may not be able to use all available features'


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

use tqtec, only: input_file, &
                 verbosity

implicit none

! Local variables
integer :: ios
character(len=512) :: input_line
character(len=8) :: exec_name
logical :: ex
double precision :: dp



if (verbosity.ge.3) then
    write(*,*) 'read_input_file: determining whether to read old or new format'
endif



! Set name of executable
#ifdef COMPILE_TQTEC
    exec_name = 'tqtec'
#elif COMPILE_TQCLIM
    exec_name = 'tqclim'
#else
    exec_name = 'NAME_ERR'
#endif




! Check to make sure input file exists
inquire(file=input_file,exist=ex)
if (.not.ex) then
    write(0,*) trim(exec_name)//': could not find input file "',trim(input_file),'"'
    call error_exit(1)
endif


! Check input file can be opened
open(unit=8,file=input_file,iostat=ios)
if (ios.ne.0) then
    write(0,*) trim(exec_name)//': something went wrong trying to open input file "',trim(input_file),'"'
    call error_exit(1)
endif


! Check the input file format
! Old version is fixed format: 1st line: model title; 2nd line: model params
! New version has format "VAR=VALUE" and allows blank lines/# commented lines
! 2nd line differentiates two formats: old format only allows numbers in this line
! whereas modern format will have characters (or be blank)
read(8,'(A)') input_line
read(8,'(A)') input_line
close(8)
read(input_line,*,iostat=ios) dp ! Check whether entry is a number
if (ios.eq.0) then
    call read_input_file_old()
else
    call read_input_file_new()
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
!          1
!    45.0000   34.0000
!          1
!     5.0000    2.0000    1.0000    0.0000   9.0000
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
                 hfvar, &
                 nthicken, &
                 thicken_dat

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
    call error_exit(1)
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


! Read surface heat flow variations
read(8,*,end=1004,iostat=ios) nhfvars
if (allocated(hfvar)) then
    deallocate(hfvar)
endif
allocate(hfvar(nhfvars,2))
do i = 1,nhfvars
    read(8,'(A)',end=1104,iostat=ios) input_line
    read(input_line,*,end=1204,iostat=ios) (hfvar(i,j),j=1,2)
enddo


! Read thickening events
read(8,*,end=1005,iostat=ios) nthicken
if (allocated(thicken_dat)) then
    deallocate(thicken_dat)
endif
allocate(thicken_dat(nthicken,5))
do i = 1,nthicken
    if (i.eq.1) then
        write(0,*) 'tqtec: reading thickening event data'
        write(0,*) 'For more precise control, use modern mode input'
        write(0,*) 'Type "tqtec -f:example" for details'
    endif
    read(8,'(A)',end=1105,iostat=ios) input_line
    read(input_line,*,end=1205,iostat=ios) (thicken_dat(i,j),j=1,5)
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
1005 if (ios.ne.0) then
    write(0,*) 'tqtec: could not find any thickening event settings in input file'
endif
ios = 0

! Errors if unable to read number of specified tectonic events
1101 if (ios.ne.0) then
    write(0,*) 'tqtec: input file only specified',i-1,' of',nburial,' burial events'
    call error_exit(1)
endif
1102 if (ios.ne.0) then
    write(0,*) 'tqtec: input file only specified',i-1,' of',nuplift,' uplift events'
    call error_exit(1)
endif
1103 if (ios.ne.0) then
    write(0,*) 'tqtec: input file only specified',i-1,' of',nthrust,' thrust events'
    call error_exit(1)
endif
1104 if (ios.ne.0) then
    write(0,*) 'tqtec: input file only specified',i-1,' of',nhfvars,' heat flow variations'
    call error_exit(1)
endif
1105 if (ios.ne.0) then
    write(0,*) 'tqtec: input file only specified',i-1,' of',nthicken,' thickening events'
    call error_exit(1)
endif

! Errors if tectonic events lines are too short
1201 if (ios.ne.0) then
    write(0,*) 'tqtec: could not parse:'
    write(0,*) '"',trim(input_line),'"'
    write(0,*) 'As: TSTART  TDURATION  BURIAL  CONDUCTIVITY'
    call error_exit(1)
endif
1202 if (ios.ne.0) then
    write(0,*) 'tqtec: could not parse:'
    write(0,*) '"',trim(input_line),'"'
    write(0,*) 'As: TSTART  TDURATION  UPLIFT'
    call error_exit(1)
endif
1203 if (ios.ne.0) then
    write(0,*) 'tqtec: could not parse:'
    write(0,*) '"',trim(input_line),'"'
    write(0,*) 'As: TSTART  1/0  BASE  DEPTH  THICKNESS'
    call error_exit(1)
endif
1204 if (ios.ne.0) then
    write(0,*) 'tqtec: could not parse:'
    write(0,*) '"',trim(input_line),'"'
    write(0,*) 'As: TSTART  HEAT_FLOW'
    call error_exit(1)
endif
1205 if (ios.ne.0) then
    write(0,*) 'tqtec: could not parse:'
    write(0,*) '"',trim(input_line),'"'
    write(0,*) 'As: TSTART  TDURATION  THICKENING  TOPCRUST  THICK0'
    call error_exit(1)
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
! 45 34
! NTHICKEN=1
! 5 2 1 0 9
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
character(len=8) :: exec_name
double precision :: max_depth
logical :: isMaxDepthDefined
logical :: isLineBlank




if (verbosity.ge.3) then
    write(*,*) 'read_input_file_new: starting'
endif



! Initialize variables
ios = 0
iend = 0
isMaxDepthDefined = .false.
isActionDefined = .false.
t_total = -1.0d99
t_geotherm_output = -1.0d99
temp_surf = -1.0d99
hf_surf = -1.0d99
cond_base = -1.0d99
hp_surf = 0.0d0
hp_dep = 0.0d0
nhorizons = 0



! Set name of executable
#ifdef COMPILE_TQTEC
    exec_name = 'tqtec'
#elif COMPILE_TQCLIM
    exec_name = 'tqclim'
#else
    exec_name = 'NAME_ERR'
#endif




! Open the input file for reading in free format
open(unit=8,file=input_file,iostat=ios)
if (ios.ne.0) then
    write(0,*) trim(exec_name)//': something went wrong trying to open input file "',trim(input_file),'"'
    call error_exit(1)
endif



! Read the file in flexible format
do while (iend.eq.0)

    ! Read the line
    read(8,'(A)',end=3451,iostat=iend) input_line


    ! Anything after "#" is a comment; ignore blank lines (function isLineBlank is below)
    if (isLineBlank(input_line)) then
        cycle
    endif


    ! All variables are in the format VAR=VALUE
    i = index(input_line,'=')
    input_line(i:i) = ' '
    read(input_line,*,iostat=ios) var, value
    if (ios.ne.0) then
        write(0,*) trim(exec_name)//': something went wrong trying to read "',trim(input_line),'"'
        call error_exit(1)
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
#ifdef COMPILE_TQCLIM
        hf_surf = hf_surf/1.0d3
#endif
    elseif (var.eq.'COND_BASE'.or.var.eq.'cond_base') then
        read(value,*) cond_base
    elseif (var.eq.'HP_SURF'.or.var.eq.'hp_surf') then
        read(value,*) hp_surf
    elseif (var.eq.'HP_DEP'.or.var.eq.'hp_dep') then
        read(value,*) hp_dep
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
    elseif (var.eq.'DT'.or.var.eq.'dt') then
        read(value,*) dt


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


    !--- Tracked Horizons ---!
    elseif (var.eq.'NHORIZONS'.or.var.eq.'nhorizons') then

        read(value,*) nhorizons

        if (nhorizons.gt.0) then

            ! Initialize horizon depth array
            if (allocated(depth)) then
                deallocate(depth)
            endif
            allocate(depth(nhorizons))
            depth = 0.0d0

            ! Skip blank or commented lines
            read(8,'(A)') input_line
            do while (isLineBlank(input_line))
                read(8,'(A)') input_line
            enddo

            read(input_line,*,iostat=ios) (depth(i),i=1,nhorizons)
                                                 ! depth

            ! Check for errors during horizon data read
            if (ios.ne.0) then
                write(0,*) trim(exec_name)//': error reading horizon parameters '
                write(0,*) 'Looking for: dep1(km) dep2(km) ... depN(km)'
                write(0,*) 'Read:        ',trim(input_line)
                call error_exit(1)
            endif
        endif


    !--- Burial Events ---!
    elseif (var.eq.'NBURIAL'.or.var.eq.'nburial') then

        read(value,*) nburial

        if (nburial.gt.0) then

            ! Initialize burial array
            if (allocated(burial_dat)) then
                deallocate(burial_dat)
            endif
            allocate(burial_dat(nburial,4))
            burial_dat = 0.0d0

            ! Read burial data
            do i = 1,nburial

                ! Skip blank or commented lines
                read(8,'(A)') input_line
                do while (isLineBlank(input_line))
                    read(8,'(A)') input_line
                enddo

                read(input_line,*,iostat=ios) (burial_dat(i,j),j=1,4)
                                                 ! start duration thickness conductivity

                ! Check for errors during burial data read
                if (ios.ne.0) then
                    write(0,*) trim(exec_name)//': error reading burial parameters for burial event',i
                    write(0,*) 'Looking for: start(Ma) duration(Ma) thick(km) cond(W/mK)'
                    write(0,*) 'Read:        ',trim(input_line)
                    call error_exit(1)
                endif
            enddo

            ! We have something to do!
            isActionDefined = .true.
        endif


    !--- Uplift/Erosion Events ---!
    elseif (var.eq.'NUPLIFT'.or.var.eq.'nuplift') then

        read(value,*) nuplift

        if (nuplift.gt.0) then

            ! Initialize uplift array
            if (allocated(uplift_dat)) then
                deallocate(uplift_dat)
            endif
            allocate(uplift_dat(nuplift,3))
            uplift_dat = 0.0d0

            ! Read uplift data
            do i = 1,nuplift

                ! Skip blank or commented lines
                read(8,'(A)') input_line
                do while (isLineBlank(input_line))
                    read(8,'(A)') input_line
                enddo

                read(input_line,*,iostat=ios) (uplift_dat(i,j),j=1,3)
                                                 ! start duration thickness

                ! Check for errors during uplift data read
                if (ios.ne.0) then
                    write(0,*) trim(exec_name)//': error reading uplift parameters for uplift event',i
                    write(0,*) 'Looking for: start(Ma) duration(Ma) thick(km)'
                    write(0,*) 'Read:        ',trim(input_line)
                    call error_exit(1)
                endif

            enddo

            ! We have something to do!
            isActionDefined = .true.
        endif


    !--- Thrusting Events ---!
    elseif (var.eq.'NTHRUST'.or.var.eq.'nthrust') then

        read(value,*) nthrust

        if (nthrust.gt.0) then

            ! Initialize thrust array
            if (allocated(thrust_dat)) then
                deallocate(thrust_dat)
            endif
            allocate(thrust_dat(nthrust,5))
            thrust_dat = 0.0d0

            ! Read thrusting data
            do i = 1,nthrust

                ! Skip blank or commented lines
                read(8,'(A)') input_line
                do while (isLineBlank(input_line))
                    read(8,'(A)') input_line
                enddo

                read(input_line,*,iostat=ios) (thrust_dat(i,j),j=1,5)
                                                 ! start upper/lower thick_init depth thick_final

                ! Check for errors during thrust data read
                if (ios.ne.0) then
                    write(0,*) trim(exec_name)//': error reading thrust parameters for thrust event',i
                    write(0,*) 'Looking for: start(Ma) 1(upper)|2(lower) thick_init(km) depth(km) thick_final(km)'
                    write(0,*) 'Read:        ',trim(input_line)
                    call error_exit(1)
                endif

            enddo

            ! We have something to do!
            isActionDefined = .true.
        endif


    !--- Surface Heat Flow Changes ---!
    elseif (var.eq.'NHFVARS'.or.var.eq.'nhfvars') then

        read(value,*) nhfvars

        if (nhfvars.gt.0) then

            ! Initialize heat flow variation array
            if (allocated(hfvar)) then
                deallocate(hfvar)
            endif
            allocate(hfvar(nhfvars,2))
            hfvar = 0.0d0

            ! Read heat flow variation data
            do i = 1,nhfvars

                ! Skip blank or commented lines
                read(8,'(A)') input_line
                do while (isLineBlank(input_line))
                    read(8,'(A)') input_line
                enddo

                read(input_line,*,iostat=ios) (hfvar(i,j),j=1,2)
                                                 ! start heat_flow

                ! Check for errors during heat flow variation data read
                if (ios.ne.0) then
                    write(0,*) trim(exec_name)//': error reading heat flow variation parameters for event',i
                    write(0,*) 'Looking for: start(Ma) heatflow(mW/m2)'
                    write(0,*) 'Read:        ',trim(input_line)
                    call error_exit(1)
                endif

            enddo

            ! We have something to do!
            isActionDefined = .true.
        endif


    !--- Bulk Thickening Events ---!
    elseif (var.eq.'NTHICKEN'.or.var.eq.'nthicken') then

        read(value,*) nthicken

        if (nthicken.gt.0) then

            ! Initialize thickening variation array
            if (allocated(thicken_dat)) then
                deallocate(thicken_dat)
            endif
            allocate(thicken_dat(nthicken,5))
            thicken_dat = 0.0d0

            do i = 1,nthicken

                ! Skip blank or commented lines
                read(8,'(A)') input_line
                do while (isLineBlank(input_line))
                    read(8,'(A)') input_line
                enddo

                read(input_line,*,iostat=ios) (thicken_dat(i,j),j=1,5)
                                                 ! start duration thickening top thick0

                ! Check for errors during thickening data read
                if (ios.ne.0) then
                    write(0,*) trim(exec_name)//': error reading thickening event parameters for thickening event',i
                    write(0,*) 'Looking for: start(Ma) duration(Ma) thicken(km) crusttop(km) thick0(km)'
                    write(0,*) 'Read:        ',trim(input_line)
                    call error_exit(1)
                endif

            enddo
            ! We have something to do!
            isActionDefined = .true.
        endif

    elseif (var.eq.'THICKENHORIZONS'.or.var.eq.'thickenhorizons'.or.var.eq.'thickenHorizons') then
        if (value.eq.'0'.or.value.eq.'N'.or.value.eq.'n'.or.value.eq.'F') then
            thickenHorizons = .false.
        elseif (value.eq.'1'.or.value.eq.'Y'.or.value.eq.'y'.or.value.eq.'T') then
            thickenHorizons = .true.
        else
            write(0,*) trim(exec_name)//': thickenHorizons must be set to T or F'
            call error_exit(1)
        endif



    !--- Surface Temperature Step Changes ---!
    elseif (var.eq.'NTEMPSTEPS') then

        read(value,*) ntempsteps

        if (ntempsteps.gt.0) then

            allocate(temp_step_dat(ntempsteps,2))

            do i = 1,ntempsteps

                ! Skip blank or commented lines
                read(8,'(A)') input_line
                do while (isLineBlank(input_line))
                    read(8,'(A)') input_line
                enddo

                read(input_line,*,iostat=ios) (temp_step_dat(i,j),j=1,2) ! start temp

                ! Check for errors during temp step change data read
                if (ios.ne.0) then
                    write(0,*) trim(exec_name)//': error reading temp step change parameter for event',i
                    write(0,*) 'Looking for: start(yr) temp(C)'
                    write(0,*) 'Read:        ',trim(input_line)
                    call error_exit(1)
                endif
            enddo

            ! We have something to do!
            isActionDefined = .true.

        endif


    !--- Surface Temperature Linear Changes ---!
    elseif (var.eq.'NTEMPRAMPS') then

        read(value,*) ntempramps

        if (ntempramps.gt.0) then

            allocate(temp_ramp_dat(ntempramps,3))

            do i = 1,ntempramps

                ! Skip blank or commented lines
                read(8,'(A)') input_line
                do while (isLineBlank(input_line))
                    read(8,'(A)') input_line
                enddo

                read(input_line,*,iostat=ios) (temp_ramp_dat(i,j),j=1,3) ! start duration temp

                ! Check for errors during temp ramp change data read
                if (ios.ne.0) then
                    write(0,*) trim(exec_name)//': error reading temp ramp change parameter for event',i
                    write(0,*) 'Looking for: start(yr) duration(yr) temp(C)'
                    write(0,*) 'Read:        ',trim(input_line)
                    call error_exit(1)
                endif
            enddo

            ! We have something to do!
            isActionDefined = .true.

        endif


    else
        write(0,*) trim(exec_name)//': no variable option named "',trim(var),'"'
        call error_exit(1)
    endif



    ! Reached the end of the file, exit
    3451 if (iend.ne.0) then
        exit
    endif
enddo




! Check that necessary variables have been defined
if (t_total.lt.0.0) then
    write(0,*) trim(exec_name)//': t_total has not been defined'
    call error_exit(1)
endif
if (temp_surf.lt.0.0) then
    write(0,*) trim(exec_name)//': temp_surf has not been defined'
    call error_exit(1)
endif
if (hf_surf.lt.0.0) then
    write(0,*) trim(exec_name)//': hf_surf has not been defined'
    call error_exit(1)
endif
if (cond_base.lt.0.0) then
    write(0,*) trim(exec_name)//': cond_base has not been defined'
    call error_exit(1)
endif
if (nhorizons.le.0) then
    write(0,*) trim(exec_name)//': no horizons have been defined'
    call error_exit(1)
endif
if (temp_surf.lt.0.0) then
    write(0,*) trim(exec_name)//': temp_surf has not been defined'
    call error_exit(1)
endif
if (.not.isActionDefined) then
    write(0,*) trim(exec_name)//': no actions have been defined'
    call error_exit(1)
endif

close(8)




if (verbosity.ge.3) then
    write(*,*) '    read_input_file_new: finished'
endif



return
end subroutine


!--------------------------------------------------------------------------------------------------!


function isLineBlank(input_line)
!----
! Check whether input_line is blank or not - includes checking whether entire line is a comment
!----
implicit none
logical :: isLineBlank
character(len=*) :: input_line
integer :: i
i = index(input_line,"#")
if (i.gt.0) then
    input_line(i:len(input_line)) = ' '
endif
if (input_line.eq.' ') then
    isLineBlank = .true.
else
    isLineBlank = .false.
endif
return
end function








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
character(len=2) :: zunit
character(len=2) :: tunit
double precision :: hf_scale


#ifdef COMPILE_TQTEC
zunit = 'km'
tunit = 'Ma'
hf_scale = 1.0d0
#elif COMPILE_TQCLIM
zunit = 'm'
tunit = 'yr'
hf_scale = 1.0d3
#endif


! Open the output file
open(unit=7,file=output_file,status='unknown')


! Write the results in the format originally specified by Kevin Furlong

! Output file name
write(7,*) trim(output_file)

! Model parameters
write(7,'(F14.3,X,A)') dz,zunit
write(7,'(F14.3,X,A)') dt,tunit
write(7,110) hp_surf
write(7,110) hp_dep
write(7,110) t_total
write(7,110) diffusivity
write(7,110) temp_factor
write(7,'(I10)') nhorizons
write(7,110) 0.0 ! II(9)
write(7,110) 0.0 ! II(10)
110 format(F14.3)

! Heat flow
do j = 2,nt_total,2
    write(7,115) hf(j)*hf_scale
enddo
115 format(E14.6)

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


subroutine print_geotherm( &
    geotherm_file,         &
    itime,                 &
    time,                  &
    total_time,            &
    nnodes,                &
    temp,                  &
    dz                     &
)
!----
! Print the geotherm at the current timestep to a file. Each geotherm in the file has the format:
!
!     > #  timestep  time_since_start  time_until_end
!     temp(1) dep(1)
!     temp(2) dep(2)
!     :
!     temp(n) dep(n)
!----


implicit none


! Arguments
character(len=*), intent(in) :: geotherm_file
integer, intent(in) :: itime
double precision, intent(in) :: time
double precision, intent(in) :: total_time
integer, intent(in) :: nnodes
double precision, intent(in) :: temp(nnodes)
double precision, intent(in) :: dz

! Local variables
integer :: i
logical :: isFileOpen


! Check whether geotherm file is open, and open it if not
inquire(file=geotherm_file,opened=isFileOpen)
if (.not.isFileOpen) then
    open(unit=12,file=geotherm_file,status='unknown')
endif

! Write header for geotherm at current timestep
write(12,'(A,I10,2E14.6)') '> #', itime, time, time-total_time

! Write geotherm at current timestep
do i = 1,nnodes
    write(12,*) temp(i),dble(i)*dz
enddo


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
write(*,2002) 'nhfvars:          ',nhfvars
if (nhfvars.gt.0) then
    write(*,'(5X,2A14)') 'start(Ma)', 'hf(mW/m2)'
    do i = 1,nhfvars
        write(*,'(5X,2F14.3)') (hfvar(i,j),j=1,2)
    enddo
endif
write(*,2002) 'nthicken:         ',nthicken
if (nthicken.gt.0) then
    write(*,'(5X,2A14,2X,3A14)') 'start(Ma)', 'duration(Ma)', 'thickening(km)', 'crust_top(km)', &
                                 'thick0(km)'
    do i = 1,nthicken
        write(*,'(5X,2F14.3,2X,3F14.3)') (thicken_dat(i,j),j=1,5)
    enddo
    write(*,*) 'Moving horizons during thickening?',thickenHorizons
endif


2001 format(5X,A18,F10.3,X,A)
2002 format(5X,A18,I10,X,A)

write(*,*)

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
!------------------------------- PROGRAM INFORMATION SUBROUTINES ----------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!






subroutine usage(str)
!----
! Print error statement (if provided) and usage, then exit with error code 1
!----

implicit none

character(len=*) :: str

! Print error statement if provided as an argument
if (str.ne.'') then
    write(0,*) trim(str)
    write(0,*)
endif

! Print usage statement
#ifdef COMPILE_TQTEC
    write(0,*) 'Usage: tqtec -i|-f INPUT_FILE -o OUTPUT_FILE [-geotherm geotherm_file] [...options...]'
#elif COMPILE_TQCLIM
    write(0,*) 'Usage: tqclim -i|-f INPUT_FILE -o OUTPUT_FILE [-geotherm geotherm_file] [...options...]'
#else
    write(0,*) 'No usage statement for this pre-processor option...exiting'
    call error_exit(1)
#endif
write(0,*)
write(0,*) 'INPUTS'
write(0,*) '-f INPUT_FILE            Input model parameter file (type "tqtec -f:example" for example)'
write(0,*) '-i[nteractive]           Interactively defined model parameters (deprecated)'
write(0,*)
write(0,*) 'OUTPUTS'
write(0,*) '-o OUTPUT_FILE           Output temperature-depth-time file for specified horizons'
write(0,*) '-geotherm GEOTHERM_FILE  Geotherms (output frequency defined in INPUT_FILE)'
write(0,*) '-timing TIMING_FILE      Timing of tectonic actions'
write(0,*)
write(0,*) 'OTHER OPTIONS'
write(0,*) '-v VERBOSITY             Verbosity level'
write(0,*) '-f:example               Print example input file and exit'
write(0,*)

! Exit with error code 1
call error_exit(1)

stop
end subroutine




!--------------------------------------------------------------------------------------------------!



subroutine print_input_file_details()
!----
! Print an example input file in modern input mode and exit
!----

implicit none

write(*,'(A)')
write(*,'(A)') '###################################################################################'
write(*,'(A)') '# The model parameters for TQTec and TQClim are input via control files.'
write(*,'(A)') '# Originally, TQTec had users interactively input parameters, but that option has'
write(*,'(A)') '# been abandoned in favor of using a control file in a flexible, commentable,'
write(*,'(A)') '# modern input format. This format sets parameters (in any order) using the format'
write(*,'(A)') '# PAR=value. Comments can be added to the file by starting a line with "#", but'
write(*,'(A)') '# note that comments cannot be inserted within parameter blocks (e.g., between'
write(*,'(A)') '# NUPLIFT= and the definition of the uplift events in the following lines).'
write(*,'(A)') '#'
write(*,'(A)') '# This text can be copied or redirected to a file, and used directly as a'
write(*,'(A)') '# comprehensive example input file with explanations for each parameter.'
write(*,'(A)') '#'
write(*,'(A)') '# NOTE: The input units are different for TQTec and TQClim:'
write(*,'(A)') '# UNITS      TQTEC       TQCLIM'
write(*,'(A)') '# Distance   kilometers  meters'
write(*,'(A)') '# Time       Ma          years'
write(*,'(A)') '###################################################################################'
write(*,'(A)')
write(*,'(A)')
write(*,'(A)') '# Example TQTec/TQClim input file in modern format'
write(*,'(A)')
write(*,'(A)')
write(*,'(A)') '#---------- MODEL TIMING PARAMETERS ----------#'
write(*,'(A)') '# Total model time (Ma|yr)'
write(*,'(A)') 'T_TOTAL=50'
write(*,'(A)')
write(*,'(A)') '# Time step size (Ma|yr)'
write(*,'(A)') 'DT=0.05'
write(*,'(A)')
write(*,'(A)') '# Time between geotherm outputs (Ma|yr)'
write(*,'(A)') 'T_GEOTHERM_OUTPUT=1   # <== Only used if -geotherm FILE is specified'
write(*,'(A)')
write(*,'(A)')
write(*,'(A)') '#---------- MODEL SPATIAL PARAMETERS ----------#'
write(*,'(A)') '# Maximum model depth (km|m)'
write(*,'(A)') 'MAX_DEPTH=100'
write(*,'(A)')
write(*,'(A)') '# Vertical node spacing (km|m)'
write(*,'(A)') 'DZ=0.1'
write(*,'(A)')
write(*,'(A)') '# Alternatively, the number of nodes can be specified and DZ will be calculated'
write(*,'(A)') 'NNODES=1000'
write(*,'(A)')
write(*,'(A)')
write(*,'(A)') '#---------- BOUNDARY CONDITIONS ----------#'
write(*,'(A)') '# By default, tqtec and tqclim have constant surface temperatures, but there are'
write(*,'(A)') '# options to have surface temperature vary over time.'
write(*,'(A)')
write(*,'(A)') '# Surface temperature (C)'
write(*,'(A)') 'TEMP_SURF=10     # <== THIS IS REQUIRED!!! (Even if you want varying surface temps)'
write(*,'(A)')
write(*,'(A)') '# Step temperature changes'
write(*,'(A)') '# Start_Time(Ma|yr) Temperature_Change(C)'
write(*,'(A)') 'NTEMPSTEPS=1'
write(*,'(A)') '0.5 0.5'
write(*,'(A)')
write(*,'(A)') '# Ramp temperature changes'
write(*,'(A)') '# Start_Time(Ma|yr) Duration(tqtec:Ma,tqclim:yr) Temperature_Change(C)'
write(*,'(A)') 'NTEMPRAMPS=1'
write(*,'(A)') '1.5 1.0 0.5'
write(*,'(A)')
write(*,'(A)')
write(*,'(A)') '# At the base of the model, TQTec and TQClim have a heat flow (not temperature!)'
write(*,'(A)') '# boundary condition. However, the boundary condition input value is the heat flow'
write(*,'(A)') '# at the surface. This is the sum of the heat flow from deep in the Earth and heat'
write(*,'(A)') '# produced by radioactive decay within the model (see HP_SURF and HP_DEP).'
write(*,'(A)')
write(*,'(A)') '# Surface heat flow (mW/m^2 = kW/km^2)'
write(*,'(A)') 'HF_SURF=75'
write(*,'(A)')
write(*,'(A)')
write(*,'(A)') '# Radioactive heat production in TQTec/TQClim is implemented as exponentially'
write(*,'(A)') '# decaying from the surface downward. This requires a surface heat production value'
write(*,'(A)') '# and a depth scale for the exponential decrease.'
write(*,'(A)')
write(*,'(A)') '# Surface heat production (uW/m^3)'
write(*,'(A)') 'HP_SURF=0'
write(*,'(A)')
write(*,'(A)') '# Heat production e-folding depth (km|m)'
write(*,'(A)') 'HP_DEP=5'
write(*,'(A)')
write(*,'(A)')
write(*,'(A)') '#---------- MATERIAL PROPERTIES ----------#'
write(*,'(A)') '# Thermal conductivity at the base of the model (W/m/K)'
write(*,'(A)') 'COND_BASE=3.0'
write(*,'(A)')
write(*,'(A)') '# Crustal layers with different conductivities'
write(*,'(A)') '# Layer_Top(km|m) Layer_Thickness(km|m) Conductivity(W/m/K)'
write(*,'(A)') 'NLAYERS=1'
write(*,'(A)') '0 10 2.5'
write(*,'(A)')

call error_exit(1)


return


end subroutine





