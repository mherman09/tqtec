!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!---------------------------------- INPUT SUBROUTINES ---------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!


subroutine usage(str)
!----
! Print error statement (if provided) and tqtec usage, then exit with error code 1
!----

implicit none

character(len=*) :: str

! Print error statement if provided as an argument
if (str.ne.'') then
    write(0,*) trim(str)
    write(0,*)
endif

! Print tqtec usage statement
write(0,*) 'Usage: tqtec -i|-f INPUT_FILE  [-o OUTPUT_FILE] [-geotherm geotherm_file] [-timing TIMING_FILE]'
write(0,*)
write(0,*) 'INPUTS'
write(0,*) '-i[nteractive]           Interactively defined model parameters'
write(0,*) '-f INPUT_FILE            Input model parameter file (type "tqtec -fd" for details)'
write(0,*)
write(0,*) 'OUTPUTS'
write(0,*) '-o OUTPUT_FILE           Output temperature-depth-time file for specified horizons'
write(0,*) '-geotherm GEOTHERM_FILE  Geotherms (output frequency defined in INPUT_FILE)'
write(0,*) '-timing TIMING_FILE      Timing of tectonic actions'
write(0,*) '-v VERBOSITY             Verbosity level'
write(0,*)
write(0,*) 'MODEL PARAMETERS'
write(0,*) 'Note: values set in INPUT_FILE new format overrides these command line values'
write(0,*) 'NNODES=<nnodes>          Number of nodes (default: 5000)'
write(0,*) 'DZ=<dz>                  Node spacing (default: 0.01 km = 10 m)'
write(0,*) 'DT=<dt>                  Time step size (default: 0.001 Ma = 1 ka)'
write(0,*)

! Exit with error code 1
call error_exit(1)

stop
end subroutine


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
character(len=512) arg
integer :: i, j, ios, narg


! Initialize control variables
ios = 0


! Initialize default values
input_file = ''
input_mode = 'user'
output_file = ''
geotherm_file = ''
timing_file = ''


! Count arguments, then exit with usage statement if no arguments
narg = command_argument_count()
if (narg.eq.0) then
    call usage('')
endif


! Parse command line arguments
i = 1
do while (i.le.narg)

    call get_command_argument(i,arg)


    if (arg.eq.'-f') then
        input_mode = 'file'
        i = i + 1
        call get_command_argument(i,input_file,status=ios)
    elseif (arg.eq.'-fd') then
        call print_input_file_details()


    elseif (arg.eq.'-i'.or.arg.eq.'-interactive') then
        input_mode = 'user'


    elseif (arg.eq.'-o') then
        i = i + 1
        call get_command_argument(i,output_file,status=ios)


    elseif (arg.eq.'-geotherm') then
        i = i + 1
        call get_command_argument(i,geotherm_file,status=ios)


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
        write(0,*) 'tqtec: could not understand response "',trim(reply),'"...'
        write(0,*) 'Exiting tqtec'
        call error_exit(1)
    endif

elseif (input_mode.eq.'file') then  ! Input file mode (for batch processing)

    call read_input_file()

else

    write(0,*) 'tqtec: no input mode named "',trim(input_mode),'"'
    call error_exit(1)

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

use tqtec, only: input_file, &
                 verbosity

implicit none

! Local variables
integer :: ios
character(len=512) :: input_line
logical :: ex
double precision :: dp



if (verbosity.ge.3) then
    write(*,*) '    read_input_file: determining whether to read old or new format'
endif



! Check to make sure input file exists
inquire(file=input_file,exist=ex)
if (.not.ex) then
    write(0,*) 'tqtec: could not find input file "',trim(input_file),'"'
    call error_exit(1)
endif


! Check input file can be opened
open(unit=8,file=input_file,iostat=ios)
if (ios.ne.0) then
    write(0,*) 'tqtec: something went wrong trying to open input file "',trim(input_file),'"'
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
        write(0,*) 'Type "tqtec -fd" for details'
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
double precision :: max_depth
logical :: isMaxDepthDefined
logical :: isLineBlank




if (verbosity.ge.3) then
    write(*,*) '    read_input_file_new: starting'
endif



! Initialize variables
ios = 0
iend = 0
isMaxDepthDefined = .false.
isTectonicActionDefined = .false.
t_total = -1.0d99
t_geotherm_output = -1.0d99
temp_surf = -1.0d99
hf_surf = -1.0d99
cond_base = -1.0d99
hp_surf = 0.0d0
hp_dep = 0.0d0
nhorizons = 0



! Open the input file for reading in free format
open(unit=8,file=input_file,iostat=ios)
if (ios.ne.0) then
    write(0,*) 'tqtec: something went wrong trying to open input file "',trim(input_file),'"'
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
        write(0,*) 'tqtec: something went wrong trying to read "',trim(input_line),'"'
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
                write(0,*) 'tqtec: error reading horizon parameters '
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
                    write(0,*) 'tqtec: error reading burial parameters for burial event',i
                    write(0,*) 'Looking for: start(Ma) duration(Ma) thick(km) cond(W/mK)'
                    write(0,*) 'Read:        ',trim(input_line)
                    call error_exit(1)
                endif
            enddo

            ! We have something to do!
            isTectonicActionDefined = .true.
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
                    write(0,*) 'tqtec: error reading uplift parameters for uplift event',i
                    write(0,*) 'Looking for: start(Ma) duration(Ma) thick(km)'
                    write(0,*) 'Read:        ',trim(input_line)
                    call error_exit(1)
                endif

            enddo

            ! We have something to do!
            isTectonicActionDefined = .true.
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
                    write(0,*) 'tqtec: error reading thrust parameters for thrust event',i
                    write(0,*) 'Looking for: start(Ma) 1(upper)|2(lower) thick_init(km) depth(km) thick_final(km)'
                    write(0,*) 'Read:        ',trim(input_line)
                    call error_exit(1)
                endif

            enddo

            ! We have something to do!
            isTectonicActionDefined = .true.
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
                    write(0,*) 'tqtec: error reading heat flow variation parameters for event',i
                    write(0,*) 'Looking for: start(Ma) heatflow(mW/m2)'
                    write(0,*) 'Read:        ',trim(input_line)
                    call error_exit(1)
                endif

            enddo

            ! We have something to do!
            isTectonicActionDefined = .true.
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
                    write(0,*) 'tqtec: error reading thickening event parameters for thickening event',i
                    write(0,*) 'Looking for: start(Ma) duration(Ma) thicken(km) crusttop(km) thick0(km)'
                    write(0,*) 'Read:        ',trim(input_line)
                    call error_exit(1)
                endif

            enddo
            ! We have something to do!
            isTectonicActionDefined = .true.
        endif

    elseif (var.eq.'THICKENHORIZONS'.or.var.eq.'thickenhorizons'.or.var.eq.'thickenHorizons') then
        if (value.eq.'0'.or.value.eq.'N'.or.value.eq.'n'.or.value.eq.'F') then
            thickenHorizons = .false.
        elseif (value.eq.'1'.or.value.eq.'Y'.or.value.eq.'y'.or.value.eq.'T') then
            thickenHorizons = .true.
        else
            write(0,*) 'tqtec: thickenHorizons must be set to T or F'
            call error_exit(1)
        endif


    else
        write(0,*) 'tqtec: no variable option named "',trim(var),'"'
        call error_exit(1)
    endif


    ! Reached the end of the file, exit
    3451 if (iend.ne.0) then
        exit
    endif
enddo


! Check that necessary variables have been defined
if (t_total.lt.0.0) then
    write(0,*) 'tqtec: t_total has not been defined'
    call error_exit(1)
endif
if (temp_surf.lt.0.0) then
    write(0,*) 'tqtec: temp_surf has not been defined'
    call error_exit(1)
endif
if (hf_surf.lt.0.0) then
    write(0,*) 'tqtec: hf_surf has not been defined'
    call error_exit(1)
endif
if (cond_base.lt.0.0) then
    write(0,*) 'tqtec: cond_base has not been defined'
    call error_exit(1)
endif
if (nhorizons.le.0) then
    write(0,*) 'tqtec: no horizons have been defined'
    call error_exit(1)
endif
if (temp_surf.lt.0.0) then
    write(0,*) 'tqtec: temp_surf has not been defined'
    call error_exit(1)
endif
if (.not.isTectonicActionDefined) then
    write(0,*) 'tqtec: no tectonic actions have been defined'
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


subroutine print_input_file_details()
implicit none
write(0,*) ''
write(0,*) 'The model parameters for TQTec are input via user interaction or control files.'
write(0,*) 'Selecting the "-i" option will lead to a series of prompts that allow you to set'
write(0,*) 'model parameters interactively. Selecting the "-f" option requires an additional'
write(0,*) 'argument defining the file name. This file must be in one of two formats:'
write(0,*)
write(0,*)
write(0,*) '1. Original TQTec Input File Format'
write(0,*) 'This format was originally fixed format, although TQTec no longer requires'
write(0,*) 'parameters to have strict widths, only that they be in the right order:'
write(0,*)
write(0,*) ' <input_file>'
write(0,*) ' <total_time_ma> <geotherm_output_ma> <temp_surf_c> <hf_surf_mW/m2> <cond_base_W/mk> <hp_surf_uW/m3> [<hp_dep_km>]'
write(0,*) ' <nlayers>'
write(0,*) ' <layer_1_dep_top_km> <layer_1_thick_km> <layer_1_cond_W/mk>'
write(0,*) ' :'
write(0,*) ' <horizon_1_dep_km> <horizon_2_dep_km> ...'
write(0,*) ' <nburial>'
write(0,*) ' <burial_1_start_ma> <burial_1_duration_ma> <burial_1_thick_km> <burial_1_cond_W/mk>'
write(0,*) ' :'
write(0,*) ' <nuplift>'
write(0,*) ' <uplift_1_start_ma> <uplift_1_duration_ma> <uplift_1_thick_km>'
write(0,*) ' :'
write(0,*) ' <nthrust>'
write(0,*) '  <thrust_1_start_ma> <upper1_or_lower2> <thrust_1_base_km> <thrust_1_depth_km> <thrust_1_thick_km>'
write(0,*) ' :'
write(0,*) ' <nthicken>'
write(0,*) '  <thicken_1_start_ma> <thicken_1_duration_ma> <thicken_1_amount_km> <thicken_1_top_km> <thicken_1_thick0_km>'
write(0,*) ' :'
write(0,*)
write(0,*)
write(0,*) '2. Modern TQTec Input File Format'
write(0,*) 'This format sets parameters (in any order) using the format PAR=value. Comments can'
write(0,*) 'be added to the file by starting a line with "#", but note that comments cannot'
write(0,*) 'be inserted within parameter blocks (e.g., between NUPLIFT= and the definition of'
write(0,*) 'the uplift events in the following lines)'
write(0,*)
write(0,*) 'T_TOTAL=total_time                      [Ma]'
write(0,*) 'T_GEOTHERM_OUTPUT=time_geotherm_output  [Ma]'
write(0,*) 'TEMP_SURF=temp_surf                     [Celsius]'
write(0,*) 'HF_SURF=hf_surf                         [mW/m2]'
write(0,*) 'COND_BASE=cond_base                     [W/mK]'
write(0,*) 'HP_SURF=hp_surf                         [uW/m3]'
write(0,*) 'HP_DEP=hp_dep                           [km]'
write(0,*) '# Define bottom of model with number of nodes or depth'
write(0,*) 'NNODES=nnodes'
write(0,*) 'MAX_DEPTH=max_depth                     [km]'
write(0,*) '# Node spacing'
write(0,*) 'DZ=dz                                   [km]'
write(0,*) 'NLAYERS=nlayers'
write(0,*) '<layer_1_dep_top_km> <layer_1_thick_km> <layer_1_cond_W/mk>'
write(0,*) ' :'
write(0,*) '# Depth horizons to track'
write(0,*) 'NHORIZONS=nhorizons'
write(0,*) '<horizon_1_dep_km> <horizon_2_dep_km> ...'
write(0,*) '# Burial events'
write(0,*) 'NBURIAL=nburial'
write(0,*) '<burial_1_start_ma> <burial_1_duration_ma> <burial_1_thick_km> <burial_1_cond_W/mk>'
write(0,*)  ':'
write(0,*) '# Uplift/erosion events'
write(0,*) 'NUPLIFT=nuplift'
write(0,*) '<uplift_1_start_ma> <uplift_1_duration_ma> <uplift_1_thick_km>'
write(0,*) ':'
write(0,*) '# Thrust faulting events (1: track in uhanging wall; 2: track in footwall)'
write(0,*) 'NTHRUST=nthrust'
write(0,*) '<thrust_1_start_ma> <upper1_or_lower2> <thrust_1_base_km> <thrust_1_depth_km> <thrust_1_thick_km>'
write(0,*) ':'
write(0,*) '# Surface heat flow variations'
write(0,*) 'NHFVARS=nhfvars'
write(0,*) '<hf_var_1_start_ma> <hf_var_1_value_mW/m2>'
write(0,*) ':'
write(0,*) ' NTHICKEN=nthicken'
write(0,*) '  <thicken_1_start_ma> <thicken_1_duration_ma> <thicken_1_amount_km> <thicken_1_top_km> <thicken_1_thick0_km>'
write(0,*) ' :'
write(0,*) 'thickenHorizons=[T|F]'
write(0,*)
write(0,*)
write(0,*) 'Example of the same model in the different formats:'
write(0,*)
write(0,*) 'Original Format:'
write(0,*) ' tqtec.in'
write(0,*) '         50         5    0.0000   30.0000    3.0000    0.0000'
write(0,*) '          0'
write(0,*) '   2.0000  4.0000  6.0000  8.0000 10.0000 12.0000 14.0000 16.0000 18.0000 20.0000'
write(0,*) '          1'
write(0,*) '    10.0000   10.0000    5.0000    2.0000'
write(0,*) '          1'
write(0,*) '    20.0000   20.0000   10.0000'
write(0,*) '          1'
write(0,*) '    40.0000         1   25.0000    0.0000  25.0000'
write(0,*) '          1'
write(0,*) '     5.0000    2.0000    1.0000    0.0000   9.0000'
write(0,*)
write(0,*)
write(0,*) 'Modern Format:'
write(0,*) 'T_TOTAL=50'
write(0,*) 'T_GEOTHERM_OUTPUT=5'
write(0,*) 'TEMP_SURF=0'
write(0,*) 'HF_SURF=30'
write(0,*) 'COND_BASE=3'
write(0,*) 'HP_SURF=0'
write(0,*) 'HP_DEP=0'
write(0,*) 'NLAYERS=0'
write(0,*) 'NHORIZONS=10'
write(0,*) '2 4 6 8 10 12 14 16 18 20'
write(0,*) 'NBURIAL=1'
write(0,*) '10 10 5 2'
write(0,*) 'NUPLIFT=2'
write(0,*) '20 20 10'
write(0,*) '45 2 1'
write(0,*) 'NTHRUST=1'
write(0,*) '40 1 25 0 25'
write(0,*) 'NHFVARS=1'
write(0,*) '15 10'
write(0,*) 'NTHICKEN=1'
write(0,*) '5 2 1 0 9   [For additional controls on thickening parameters, see Modern Format above]'




write(0,*) ''
call error_exit(1)
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





