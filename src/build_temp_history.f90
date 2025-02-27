!--------------------------------------------------------------------------------------------------!
! build_temp_history                                                                               !
!                                                                                                  !
! Authors:                                                                                         !
!     - Matt Herman                                                                                !
!                                                                                                  !
! build_temp_history allows the user to construct a piecewise linear temperature history file in   !
! the format output by readtqtec, e.g., for input to minage.                                       !
!--------------------------------------------------------------------------------------------------!

module build_temp_history

    ! Output file
    character(len=512) :: output_file

    ! Input file
    character(len=512) :: input_temp_file

    ! Variables for constructing temperature history
    double precision :: total_time     ! Total time (Ma)
    double precision :: dt             ! Timestep (Ma)
    integer :: npts                    ! Number of points in timeseries
    double precision :: time0(1000)    ! Times to set temperatures
    double precision :: temp0(1000)    ! Temperatures at defined times

end module

!==================================================================================================!

program main
!----
! Manually construct temperature file in readtqtec output format
!----

use build_temp_history

implicit none

! Local variables
integer :: i
integer :: j
integer :: n
integer :: ios
double precision :: t
double precision :: temp
character(len=512) :: input_line

! Parse command line
call gcmdln()

! Calculate number of points in time series
npts = int(total_time/dt)

! Set initial temperature at time=0
time0(1) = 0.0d0
temp0(1) = 0.0d0


! Read temperatures from file or stdin
if (input_temp_file.ne.'') then

    ! Open temperature file and read times/temperatures
    open(unit=11,file=input_temp_file,status='old')
    n = 0
    do
        n = n + 1
        read(11,*,end=101,iostat=ios) time0(n),temp0(n)
        101 if (ios.ne.0) then
            n = n - 1
            exit
        endif
    enddo
    close(11)

    ! Make sure that there is an initial temperature
    if (time0(1).gt.0.0d0) then
        write(0,*) 'build_temp_history [WARNING]: temp not set at time=0 Ma in file "'&
                   //trim(input_temp_file)//'"'
        write(0,*) 'build_temp_history [WARNING]: setting initial temp=0 C'
        do i = n,1,-1
            time0(i+1) = time0(i)
            temp0(i+1) = temp0(i)
        enddo
        n = n + 1
        time0(1) = 0.0d0
        temp0(1) = 0.0d0
    endif

else

    ! Read times/temperatures from standard input
    time0(1) = 0.0d0
    write(*,*) 'Starting Temp(C)?'
    read(*,*) temp0(1)
    n = 1
    do
        write(*,*) 'Time(Ma) Temp(C)? [Type "end" to exit]'
        read(*,'(A)') input_line
        if (input_line.eq.'end') then
            exit
        endif
        n = n + 1
        read(input_line,*) time0(n),temp0(n)
    enddo
endif

! Set final time0 to be the ending time (only used if input time is less than final model time)
time0(n+1) = total_time

! Sanity checks
if (abs(time0(1)-0.0d0).gt.1.0d-8) then
    write(0,*) 'build_temp_history [ERROR]: starting time is somehow not 0 Ma'
    stop
endif
do i = 2,n
    if (time0(i).le.time0(i-1)) then
        write(0,*) 'build_temp_history [ERROR]: input times must increase monotonically'
        stop
    endif
enddo

! Write fixed time-temp coordinates
write(*,'(A)') 'Generating piecewise linear temperature history...'
write(*,'(A)') ' Time0(Ma)  Temp0(C)'
do i = 1,n
    write(0,'(F10.3,F10.3)') time0(i),temp0(i)
enddo

! Write output file
open(file=output_file,unit=12,status='unknown')
write(12,'(2F10.3,I6)') -total_time, dt, npts
i = 1
do j = 1,npts
    t = dble(j)*dt
    if (t.gt.time0(i+1)) then
        i = i + 1
    endif
    if (i.eq.n) then
        temp = temp0(i)
    else
        temp = temp0(i) + (t-time0(i))*(temp0(i+1)-temp0(i))/(time0(i+1)-time0(i))
    endif
    write(12,'(1PE14.6)') temp
enddo
close(12)
end

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()
use build_temp_history
implicit none
! Local variables
character(len=512) :: arg
integer :: i, narg
total_time = 0.0d0
dt = 0.0d0
output_file = ''
input_temp_file = ''
! Count arguments, then exit with usage statement if no arguments
narg = command_argument_count()
if (narg.eq.0) then
    call usage('')
endif
! Parse command line arguments
i = 1
do while (i.le.narg)
    call get_command_argument(i,arg)
    if (arg.eq.'-o') then
        i = i + 1
        call get_command_argument(i,output_file)
    elseif (arg.eq.'-total_time') then
        i = i + 1
        call get_command_argument(i,arg)
        read(arg,*) total_time
    elseif (arg.eq.'-dt') then
        i = i + 1
        call get_command_argument(i,arg)
        read(arg,*) dt
    elseif (arg.eq.'-input_temp_file') then
        i = i + 1
        call get_command_argument(i,input_temp_file)
    else
        call usage('build_temp_history [ERROR]: no option "'//trim(arg)//'"')
    endif
    i = i + 1
enddo
if (output_file.eq.'') then
    call usage('build_temp_history [ERROR]: set output file with -o')
endif
if (total_time.le.0.0d0) then
    call usage('build_temp_history [ERROR]: set total time with -total_time')
endif
if (dt.le.0.0d0) then
    call usage('build_temp_history [ERROR]: set dt with -dt')
endif
return
end

!--------------------------------------------------------------------------------------------------!

subroutine usage(str)
implicit none
character(len=*) :: str
! Print error statement if provided as an argument
if (str.ne.'') then
    write(0,*) trim(str)
    write(0,*)
endif
write(0,*) 'Usage: build_temp_history -o FILE -total_time T_MA -dt DT_MA [-input_temp_file FILE]'
stop
end subroutine
