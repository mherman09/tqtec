module build_temp_history
    double precision :: total_time
    double precision :: dt
    integer :: npts
    double precision :: time0(1000)
    double precision :: temp0(1000)
    character(len=512) :: temp_file
end module

program main
use build_temp_history
implicit none
integer :: i, j, n, ios
double precision :: t, temp
character(len=512) :: input_line
call gcmdln()
npts = int(total_time/dt)
time0(1) = 0.0d0
temp0(1) = 0.0d0
if (temp_file.ne.'') then
    open(unit=11,file=temp_file,status='old')
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
else
    time0(1) = 0.0d0
    write(0,*) 'Starting Temp(C)?'
    read(*,*) temp0(1)
    n = 1
    do
        write(0,*) 'Time(Ma) Temp(C)? [Type "end" to exit]'
        read(*,'(A)') input_line
        if (input_line.eq.'end') then
            exit
        endif
        n = n + 1
        read(input_line,*) time0(n),temp0(n)
    enddo
endif
time0(n+1) = total_time

! WRITE RESULTS
write(*,'(2F10.3,I6)') -total_time, dt, npts
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
    write(*,'(1PE14.6)') temp
enddo
end

subroutine gcmdln()
use build_temp_history
implicit none
! Local variables
character(len=512) :: arg
integer :: i, narg
total_time = 0.0d0
dt = 0.0d0
temp_file = ''
! Count arguments, then exit with usage statement if no arguments
narg = command_argument_count()
if (narg.eq.0) then
    call usage('')
endif
! Parse command line arguments
i = 1
do while (i.le.narg)
    call get_command_argument(i,arg)
    if (arg.eq.'-total_time') then
        i = i + 1
        call get_command_argument(i,arg)
        read(arg,*) total_time
    elseif (arg.eq.'-dt') then
        i = i + 1
        call get_command_argument(i,arg)
        read(arg,*) dt
    elseif (arg.eq.'-temp_file') then
        i = i + 1
        call get_command_argument(i,temp_file)
    endif
    i = i + 1
enddo
if (total_time.le.0.0d0) then
    write(0,*) 'build_temp_history [ERROR]: set total time with -total_time'
    call usage('')
endif
if (dt.le.0.0d0) then
    write(0,*) 'build_temp_history [ERROR]: set dt with -dt'
    call usage('')
endif
return
end

subroutine usage(str)
implicit none
character(len=*) :: str
! Print error statement if provided as an argument
if (str.ne.'') then
    write(0,*) trim(str)
    write(0,*)
endif
write(0,*) 'Usage: build_temp_history -total_time T_MA -dt DT_MA [-temp_file FILE]'
stop
end subroutine
