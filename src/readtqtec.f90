!----
! readTQTec
!
! Authors:
!     - Kevin Furlong (original Fortran 77 program)
!     - Matt Herman (Modern Fortran version, i.e., what you are looking at right now!)
!     - Chris Guzofski, Matt Legg, Rachel Piotraschke (PSU MS theses)
!
! C     READS THE OUTPUT FROM THE PROGRAM TQTEC
!----


module readtqtec


! Input file
character(len=512) :: tqtec_output_file            ! tqtec output with temp/dep/time data

! Outputs
character(len=512) :: temp_file                    ! temperature-timestep file
character(len=512) :: dep_file                     ! depth-timestep file
character(len=512) :: time_file                    ! time-timestep file
character(len=512) :: hf_file                      ! heat flow-timestep file
character(len=512) :: tti_file                     ! tti-timestep file
character(len=512) :: vr_file                      ! vitrinite reflectance-timestep file
character(len=512) :: closure_file                 ! closure temp-depth file
integer :: nclosure                                ! number of closure temperatures to track
double precision, allocatable :: closure_temps(:)  ! closure temperatures to track
logical, allocatable :: isPartialAnnealTemp(:)     ! flags for specifying partial annealing zone
character(len=512) :: temp_range_file              ! readtqtec output temperature range file
integer :: ntemprange                              ! number of temperature ranges to track
double precision, allocatable :: temp_ranges(:,:)  ! temperature ranges


end module


!==================================================================================================!


program main
!----
! Reads raw output from TQTec and prints results for post-processing
!----


use readtqtec


implicit none

! Parameters from TQTec model
integer :: nt_total                               ! number of time steps                            Q1 (updated), II(5)
double precision :: dz                            ! node spacing (km)                               H1, II(1)
double precision :: dt                            ! time step interval (Ma)                         K1, II(2)
double precision :: t_total                       ! total model time (Ma)                           Q1 (initial), II(5)
double precision, allocatable :: hf(:)            ! heat flow                                       Q
integer :: nhorizons                              ! number of horizons                              10
integer, allocatable :: depth_node(:)             ! horizon nodes                                   Y (updated)
double precision :: diffusivity                   ! diffusivity                                     D1, II(6)
double precision :: temp_surf                     ! surface temperature (C)                         W(1)
double precision :: hp_surf                       ! surface heat production                         A1, II(3)
double precision :: hp_dep                        ! depth of heat production                        B1, II(4)

! TQTec output array
double precision, allocatable :: results(:,:,:)   ! temperature and depth for each timestep/depth   r1

! Local variables
integer :: i
integer :: j
integer :: k
integer :: l
integer :: ios
character(len=1) :: dummy
character(len=32) :: fmt_string
logical :: ex, isClosed, isInRange
double precision :: xmin, closure_time, final_depth, temp1
double precision, allocatable :: tti_exp(:)
double precision, allocatable :: tti(:)
double precision, allocatable :: vr0(:)




! Parse command line options
call gcmdln()



! Check whether tqtec output file exists
! NOTE: subroutine output() in tqtec_io.f90 produces this file
inquire(file=tqtec_output_file,exist=ex)
if (.not.ex) then
    call usage('readtqtec: no file found named "'//trim(tqtec_output_file)//'"')
endif



! Open the file and start reading
open(unit=8,file=tqtec_output_file,status='old')



! First line is just the file name
read(8,'(A)') dummy


! Model parameters
read(8,*,iostat=ios,err=1001,end=2001) dz
read(8,*,iostat=ios,err=1001,end=2001) dt
read(8,*,iostat=ios,err=1001,end=2001) hp_surf
read(8,*,iostat=ios,err=1001,end=2001) hp_dep
read(8,*,iostat=ios,err=1001,end=2001) t_total
read(8,*,iostat=ios,err=1001,end=2001) diffusivity
read(8,*,iostat=ios,err=1001,end=2001) temp_surf ! temp_factor?
read(8,*,iostat=ios,err=1001,end=2001) nhorizons
read(8,'(A)') dummy
read(8,'(A)') dummy

! Check for read errors or end-of-file
1001 if (ios.ne.0) then
    write(0,*) 'readtqtec: error reading variable in model parameter block'
    call error_exit(1)
endif
2001 if (ios.ne.0) then
    write(0,*) 'readtqtec: reached end of file before end of model parameter block'
    call error_exit(1)
endif


! Set starting time to be (negative) model run time (Ma)
xmin = t_total

! Total number of timesteps
nt_total = int(t_total/(2.0d0*dt))



! Read heat flow at each timestep
if (allocated(hf)) then
    deallocate(hf)
endif
allocate(hf(nt_total))

do j = 1,nt_total
    read(8,*,iostat=ios,err=1002,end=2002) hf(j)
enddo

! Check for read errors or end-of-file
1002 if (ios.ne.0) then
    write(0,*) 'readtqtec: error reading heat flow'
    call error_exit(1)
endif
2002 if (ios.ne.0) then
    write(0,*) 'readtqtec: reached end of file before end of heat flow values'
    call error_exit(1)
endif



! For each tracked horizon, read temperature and depth at each timestep
if(allocated(results)) then
    deallocate(results)
endif
allocate(results(nt_total,2,nhorizons))

do k = 1,nhorizons
    do j = 1,2
        do i = 1,nt_total
            read(8,*,err=1003,end=2003) results(i,j,k)
        enddo
    enddo
enddo

! Check for read errors or end-of-file
1003 if (ios.ne.0) then
    write(0,*) 'readtqtec: error reading temperature and depth of tracked horizons'
    call error_exit(1)
endif
2003 if (ios.ne.0) then
    write(0,*) 'readtqtec: reached end of file before end of horizons temperatures and depths'
    call error_exit(1)
endif



! Read final depth of each horizon
if (allocated(depth_node)) then
    deallocate(depth_node)
endif
allocate(depth_node(nhorizons))

do i = 1,nhorizons
    read(8,*,err=1004,end=2004) depth_node(i)
enddo

! Check for read errors or end-of-file
1004 if (ios.ne.0) then
    write(0,*) 'readtqtec: error reading horizon depths'
    call error_exit(1)
endif
2004 if (ios.ne.0) then
    write(0,*) 'readtqtec: reached end of file before end of horizon depths'
    call error_exit(1)
endif




! Close the TQTec output file
close(8)








!**************************************************************************************************!
!*** PRINT REQUESTED OUTPUTS **********************************************************************!
!**************************************************************************************************!


! Temperature history
if (temp_file.ne.'') then

    open(unit=9,file=temp_file,status='unknown')

    ! Header: time_before_present(Ma)  sample_rate(Ma)  number_of_points
    write(9,'(2F10.3,I6)') -xmin, 2.0d0*dt, nt_total
            ! sample_rate = 2*dt because every 2 pts are saved in tqtec_io.f90 subroutine output()

    ! Write temperature for each horizon at each timestep
    write(fmt_string,'("("I6,"F8.3)")') nhorizons
    do l = 1,nt_total
        write(9,fmt_string) (results(l,1,i),i=1,nhorizons)
    enddo

    close(9)
endif



! Depth history
if (dep_file.ne.'') then

    open(unit=9,file=dep_file,status='unknown')

    ! Header: time_before_present(Ma)  sample_rate(Ma)  number_of_points
    write(9,'(2F10.3,I6)') -xmin,2.0d0*dt,nt_total 
            ! sample_rate = 2*dt because every 2 pts are saved in tqtec_io.f90 subroutine output()

    ! Write depth for each horizon at each timestep
    write(fmt_string,'("("I6,"F8.3)")') nhorizons
    do l = 1,nt_total
        write(9,fmt_string) (-1.0d0*results(l,2,i)*dz,i=1,nhorizons)
    enddo

    close(9)
endif



! Surface heat flow history
if (hf_file.ne.'') then

    open(unit=9,file=hf_file,status='unknown')

    ! Header: time_before_present(Ma)  sample_rate(Ma)  number_of_points
    write(9,'(2F10.3,I6)') -xmin,2.0d0*dt,nt_total
            ! sample_rate = 2*dt because every 2 pts are saved in tqtec_io.f90 subroutine output()

    ! Write surface heat flow at each timestep
    do i = 1,nt_total
        write(9,*) hf(i)
    enddo

    close(9)
endif



! Time history
if (time_file.ne.'') then

    open(unit=9,file=time_file,status='unknown')

    ! Header: time_before_present(Ma)  sample_rate(Ma)  number_of_points
    write(9,'(2F10.3,I6)') -xmin,2.0d0*dt,nt_total
            ! sample_rate = 2*dt because every 2 pts are saved in tqtec_io.f90 subroutine output()

    ! Write time in Ma for each horizon at each timestep
    write(fmt_string,'("("I6,"F8.3)")') nhorizons
    do l = 1,nt_total
        write(9,fmt_string) (l*dt*2.0d0,i=1,nhorizons)
    enddo

    close(9)
endif




! Time-Temperature Index of Maturity (TTI) history
if (tti_file.ne.'') then

    if (allocated(tti)) then
        deallocate(tti)
    endif
    allocate(tti(nhorizons))
    if (allocated(tti_exp)) then
        deallocate(tti_exp)
    endif
    allocate(tti_exp(nhorizons))

    open(unit=9,file=tti_file,status='unknown')

    ! Header: time_before_present(Ma)  sample_rate(Ma)  number_of_points
    write(9,'(2F10.3,I6)') -xmin,2.0d0*dt,nt_total
            ! sample_rate = 2*dt because every 2 pts are saved in tqtec_io.f90 subroutine output()

    ! Calculate and print TTI for each horizon at each timestep
    write(fmt_string,'("("I6,"E14.6)")') nhorizons
    do l = 1,nt_total
        tti_exp = 0.1d0*(results(l,1,:)-105.0d0)
        if (l.eq.1) then
            tti = 2.0d0*dt * 2**tti_exp
        else
            tti = tti + 2.0d0*dt * 2**tti_exp
        endif
        write(9,fmt_string) (tti(i),i=1,nhorizons)
    enddo

    close(9)
endif



! Vitrinite reflectance history
if (vr_file.ne.'') then

    if (allocated(tti)) then
        deallocate(tti)
    endif
    allocate(tti(nhorizons))
    if (allocated(tti_exp)) then
        deallocate(tti_exp)
    endif
    allocate(tti_exp(nhorizons))
    if (allocated(vr0)) then
        deallocate(vr0)
    endif
    allocate(vr0(nhorizons))

    open(unit=9,file=vr_file,status='unknown')

    ! Header: time_before_present(Ma)  sample_rate(Ma)  number_of_points
    write(9,'(2F10.3,I6)') -xmin,2.0d0*dt,nt_total
            ! sample_rate = 2*dt because every 2 pts are saved in tqtec_io.f90 subroutine output()

    ! Calculate TTI and VR0 and print VR0 for each horizon at each timestep
    write(fmt_string,'("("I6,"E14.6)")') nhorizons
    do l = 1,nt_total
        tti_exp = 0.1d0*(results(l,1,:)-105.0d0)
        if (l.eq.1) then
            tti = 2.0d0*dt * 2**tti_exp
        else
            tti = tti + 2.0d0*dt * 2**tti_exp
        endif
        vr0 = (tti/41.0d0) ** 0.22d0
        write(9,fmt_string) (vr0(i),i=1,nhorizons)
    enddo

    close(9)
endif


! C
! C     WRITE PRODUCTION TYPE I PLOT DATA
! C
! 540   WRITE(*,*)'TYPE NAME OF OUTPUT FILE FOR PROD TYPE I PLOTS'
!       READ(*,100) OUTFILE
!       OPEN(UNIT=9,FILE=OUTFILE)
!          WRITE(9,175) -XMIN,II(2)*2,Q1,0.0,1000.0
!       DATA E/48,50,52,54,56,58,60/,
!      *   C/8,10,20,26,810,11,13/,
!      *   A/7*3E29/
!       DO 80 L=1,5
!          DO 81 I=1,7
!             S(I)=0.0
! 81       CONTINUE
!          DO 88 ILOOP=1,Q1
!             TEMP=0
!             DO 89 M=1,7
!                S(M)=S(M)+EXP(-E(M)*1000/(R1*(R(ILOOP,1,L)+273)))
!      *            *A(M)*II(2)*2
!                TEMP=TEMP+C(M)*(1-EXP(-S(M)))
! 89          CONTINUE
!             PRODN(ILOOP)=TEMP
! 88       CONTINUE
!          DIFF(1)=0
!          DO 84 I=2,Q1
!            DIFF(I)=(PRODN(I)-PRODN(I-1))/II(2)*2
! 84       CONTINUE
!             DO 85 K=10,Q1,10
! 85             WRITE(9,160)(PRODN(I),I=K-9,K)
!             IF(K.LT.Q1) WRITE(9,160)(PRODN(I),I=K+1,Q1)
!          IF (K.NE.Q1.AND.K.NE.(Q1+10)) THEN
!             KNEW=ABS(Q1-K)
!             WRITE(9,160)(PRODN(I),I=Q1-KNEW+1,Q1)
!          ENDIF
! 80    CONTINUE
!       CLOSE(9)
!       GOTO 500
! C
! C     WRITE PRODUCTION TYPE II PLOT DATA
! C
! 550   WRITE(*,*) 'TYPE NAME OF OUTPUT FILE FOR PROD TYPE II PLOTS'
!       READ(*,100) OUTFILE
!       OPEN(UNIT=9, FILE =OUTFILE)
!       WRITE(9,175) -XMIN,II(2)*2,Q1,0.0,700.0
!       DATA E2/40,46,48,50,52,54,56,58,60/,
!      *C2/6,4,9,32,132,302,104,35,6/,
!      *A2/9*3E29/
!       DO 280 L=1,5
!          DO 281 I=1,9
!             S2(I)=0.0
! 281      CONTINUE
!          DO 288 ILOOP=1,Q1
!          TEMP=0
!          DO 289 M=1,9
!             S2(M)=S2(M)+EXP(-E2(M)*1000/(R1*(R(ILOOP,1,L)+273)))
!      *           *A2(M)*II(2)*2
!             TEMP=TEMP+C2(M)*(1-EXP(-S2(M)))
! 289      CONTINUE
!          PRODN(ILOOP)=TEMP
! 288   CONTINUE
!       DIFF(1)=0
!       DO 284 I=2,Q1
!          DIFF(I)=(PRODN(I)-PRODN(I-1))/II(2)*2
! 284   CONTINUE
!          DO 285 K=10,Q1,10
! 285      WRITE(9,160)(PRODN(I),I=K-9,K)
!          IF (K.NE.Q1.AND.K.NE.(Q1+10)) THEN
!             KNEW=ABS(Q1-K)
!             WRITE(9,160)(PRODN(I),I=Q1-KNEW+1,Q1)
!          ENDIF
! 280   CONTINUE
!       CLOSE(9)
!       GOTO 500
! C
! C     WRITE PRODUCTION TYPE III PLOT DATA
! C
! 560   WRITE(*,*) 'TYPE NAME OF OUTPUT FILE FOR PROD TYPE III PLOTS'
!       READ(*,100) OUTFILE
!       OPEN (UNIT=9,FILE=OUTFILE)
!          WRITE(9,175) -XMIN,II(2)*2,Q1,0.0,250.0
!       DATA E3/50,52,54,56,58,60,62,64,66,68,70,72,74/,
!      *C3/1,5,6,42,82,60,23,12,7,5,3,2,2/,
!      *A3/13*3E29/
!       DO 380 L=1,5
!          DO 381 I=1,13
!             S3(I)=0.0
! 381      CONTINUE
!          DO 388 ILOOP=1,Q1
!             TEMP=0
!             DO 389 M=1,13
!                S3(M)=S3(M)+EXP(-E3(M)*1000/(R1*(R(ILOOP,1,L)+273)))
!      *               *A3(M)*II(2)*2
!                TEMP=TEMP+C3(M)*(1-EXP(-S3(M)))
! 389         CONTINUE
!             PRODN(ILOOP)=TEMP
! 388      CONTINUE
!          DIFF(1)=0
!          DO 384 I=2,Q1
!             DIFF(I)=(PRODN(I)-PRODN(I-1))/II(2)*2
! 384      CONTINUE
!          DO 385 K=10,Q1,10
! 385      WRITE(9,160)(PRODN(I),I=K-9,K)
!          IF (K.NE.Q1.AND.K.NE.(Q1+10)) THEN
!             KNEW=ABS(Q1-K)
!             WRITE(9,160)(PRODN(I),I=Q1-KNEW+1,Q1)
!          ENDIF
! 380   CONTINUE
!       CLOSE(9)
!       GOTO 500



! Closure temperature depth and time
if (closure_file.ne.'') then
    open(unit=13,file=closure_file,status='unknown')
    do i = 1,nclosure

        ! Write closure temperature to header
        if (.not.isPartialAnnealTemp(i)) then
            write(13,'(A,F10.3)') '>',closure_temps(i)
        else
            write(13,'(A,F10.3,A)') '>',closure_temps(i),' # Partial annealing zone'
        endif

        ! Find time that each horizon passed through closure temperature
        do j = 1,nhorizons

            ! Initialize horizon as not closed at time zero
            isClosed = .false.

            do k = 1,nt_total
                if (results(k,1,j).lt.closure_temps(i)) then
                    ! Horizon temperature < closure temperature => CLOSED!
                    if (.not.isClosed) then
                        ! Save open-to-closed time
                        closure_time = k*dt*2.0d0
                        isClosed = .true.
                    endif
                else
                    ! Horizon temperature >= closure temperature => OPEN!
                    isClosed = .false.
                endif
            enddo

            if (isClosed) then
                final_depth = -1.0d0*results(nt_total,2,j)*dz
                write(13,*) closure_time,final_depth
            else
                write(13,'(A,I5,A)') '# horizon',j,' not closed at end of model run'
            endif
        enddo

        if (isPartialAnnealTemp(i)) then
            ! Find time that each horizon passed through previous closure temperature
            do j = nhorizons,1,-1

                ! Initialize horizon as not closed at time zero
                isClosed = .false.

                do k = 1,nt_total
                    if (results(k,1,j).lt.closure_temps(i-1)) then
                        ! Horizon temperature < closure temperature => CLOSED!
                        if (.not.isClosed) then
                            ! Save open-to-closed time
                            closure_time = k*dt*2.0d0
                            isClosed = .true.
                        endif
                    else
                        ! Horizon temperature >= closure temperature => OPEN!
                        isClosed = .false.
                    endif
                enddo

                if (isClosed) then
                    final_depth = -1.0d0*results(nt_total,2,j)*dz
                    write(13,*) closure_time,final_depth,'# previous closure temp'
                else
                    write(13,'(A,I5,A)') '# horizon',j,' not closed at end of model run'
                endif
            enddo
        endif

    enddo
    close(13)
endif


!----
! Write when each horizon is within temperature range
!----
if (temp_range_file.ne.'') then
    open(unit=14,file=temp_range_file,status='unknown')
    do i = 1,ntemprange

        ! Make sure temp1<temp2
        if (temp_ranges(i,2).lt.temp_ranges(i,1)) then
            temp1 = temp_ranges(i,1)
            temp_ranges(i,1) = temp_ranges(i,2)
            temp_ranges(i,2) = temp1
        endif

        ! Write temperature range to header
        write(14,'(A,2F10.3)') '>',temp_ranges(i,1),temp_ranges(i,2)

        ! Find time that each horizon was in temperature range
        do j = 1,nhorizons

            ! Initialize horizon as not in range at time zero
            isInRange = .false.

            final_depth = -1.0d0*results(nt_total,2,j)*dz

            do k = 1,nt_total
                temp1 = results(k,1,j)
                if (temp_ranges(i,1).le.temp1 .and. temp1.le.temp_ranges(i,2)) then
                    if (.not.isInRange) then
                        ! Temperature of horizon entered range
                        closure_time = k*dt*2.0d0
                        isInRange = .true.
                        ! write(14,'(A,I6,A,F10.3)') '> Horizon',j,' temp=',temp1
                        write(14,'(A)') '>'
                        write(14,*) closure_time,final_depth
                    endif
                else
                    if (isInRange) then
                        ! Temperature of horizon exited range
                        closure_time = k*dt*2.0d0
                        isInRange = .false.
                        write(14,*) closure_time,final_depth
                    endif
                endif
            enddo

            if (isInRange) then
                write(14,*) closure_time,final_depth
            else
                write(14,'(A,I5,A)') '# horizon',j,' not in temp range at end of model run'
            endif
        enddo

    enddo
    close(14)
endif



end

!--------------------------------------------------------------------------------------------------!


subroutine gcmdln()

use readtqtec, only: tqtec_output_file, &
                     temp_file, &
                     dep_file, &
                     time_file, &
                     hf_file, &
                     tti_file, &
                     vr_file, &
                     closure_file, &
                     nclosure, &
                     closure_temps, &
                     isPartialAnnealTemp, &
                     temp_range_file, &
                     ntemprange, &
                     temp_ranges

implicit none

! Local variables
character(len=512) arg
integer :: i, j, k, ios, narg
logical :: isNumber
double precision :: dp


! Initialize control variables
ios = 0

! Initialize defaults
tqtec_output_file = ''
temp_file = ''
dep_file = ''
time_file = ''
hf_file = ''
tti_file = ''
vr_file = ''
closure_file = ''
temp_range_file = ''
nclosure = 0
isNumber = .false.


narg = command_argument_count()
if (narg.eq.0) then
    call usage('')
endif

i = 1
if (i.le.narg) then
    call get_command_argument(i,tqtec_output_file)
    if (tqtec_output_file.eq.'-format') then
        print *,'here'
        call usage('Output file formats:fileformats')
    endif
endif
i = i + 1

do while (i.le.narg)

    call get_command_argument(i,arg)


    if (arg.eq.'-temp') then
        i = i + 1
        call get_command_argument(i,temp_file)
    elseif (arg.eq.'-dep') then
        i = i + 1
        call get_command_argument(i,dep_file)
    elseif (arg.eq.'-time') then
        i = i + 1
        call get_command_argument(i,time_file)
    elseif (arg.eq.'-hf') then
        i = i + 1
        call get_command_argument(i,hf_file)
    elseif (arg.eq.'-tti') then
        i = i + 1
        call get_command_argument(i,tti_file)
    elseif (arg.eq.'-vr') then
        i = i + 1
        call get_command_argument(i,vr_file)


    elseif (arg.eq.'-closure') then
        i = i + 1
        call get_command_argument(i,closure_file)                               ! Read output file name
        do                                                                      ! Count number of closure Ts
            i = i + 1
            call get_command_argument(i,arg)
            k = index(arg,'p')                                                  ! Ignore partial annealing details
            if (k.ne.0) then
                arg(k:k) = ''
            endif
            read(arg,*,iostat=ios) dp                                           ! Read as double precision value
            if (ios.ne.0) then                                                  ! New or last arg => IOS=0
                exit                                                            ! End count
            else
                nclosure = nclosure + 1                                         ! Found closure T, add to count
            endif
        enddo
        i = i - nclosure - 1                                                    ! Reset arg index to 1st closure T
        allocate(closure_temps(nclosure))                                       ! Allocate closure T array
        allocate(isPartialAnnealTemp(nclosure))                                 ! Allocate partial annealing flag array
        closure_temps = 0.0d0
        isPartialAnnealTemp = .false.
        do j = 1,nclosure                                                       ! Read values and fill arrays
            i = i + 1
            call get_command_argument(i,arg)
            k = index(arg,'p')                                                  ! Check whether "p" flag is set
            if (k.ne.0) then                                                    ! Found partial annealing T
                isPartialAnnealTemp(j) = .true.
                arg(k:k) = ''
            endif
            read(arg,*) closure_temps(j)
        enddo

    elseif (arg.eq.'-temp:range') then
        i = i + 1
        call get_command_argument(i,temp_range_file)                            ! Read output file name
        do                                                                      ! Count number of temp ranges
            i = i + 1
            if (i.gt.narg) then
                exit
            endif
            call get_command_argument(i,arg)
            k = index(arg,'-')
            arg(k:k) = ' '
            read(arg,*,iostat=ios) dp                                           ! Read as double precision value
            if (ios.ne.0) then                                                  ! New or last arg => IOS=0
                exit                                                            ! End count
            else
                ntemprange = ntemprange + 1                                     ! Found temp range, add to count
            endif
        enddo
        i = i - ntemprange - 1                                                  ! Reset arg index to 1st temp range
        allocate(temp_ranges(ntemprange,2))                                     ! Allocate temp range array
        temp_ranges = 0.0d0
        do j = 1,ntemprange                                                     ! Read values and fill arrays
            i = i + 1
            call get_command_argument(i,arg)
            k = index(arg,'-')
            arg(k:k) = ' '
            read(arg,*) temp_ranges(j,1),temp_ranges(j,2)
        enddo
    elseif (arg.eq.'-format') then
        print *,'here'
        call usage('Output file formats:fileformats')
    endif
    i = i + 1
enddo


return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine usage(str)
implicit none
character(len=*) :: str
integer :: i
logical :: printFileFormats
printFileFormats = .false.
i = index(str,':fileformats')
if (i.ne.0) then
    printFileFormats = .true.
else
    i = len_trim(str) + 1
endif
if (str.ne.'') then
    write(0,*) trim(str(1:i-1))
    write(0,*)
endif
write(0,*) 'Usage: readtqtec TQTEC_OUTPUT_FILE'
write(0,*) '                 [-temp TEMP_FILE]'
write(0,*) '                 [-dep DEP_FILE]'
write(0,*) '                 [-time TIME_FILE]'
write(0,*) '                 [-hf HF_FILE]'
write(0,*) '                 [-tti TTI_FILE]'
write(0,*) '                 [-closure FILE  T1[p] T2[p] T3[p]...]'
write(0,*) '                 [-temp:range FILE  T1-T2]'
write(0,*) '                 [-format]'
write(0,*)
write(0,*) 'TQTEC_OUTPUT_FILE                   TQTec output file'
write(0,*)
write(0,*) '-temp TEMP_FILE                      Temperature over time for each horizon'
write(0,*) '-dep DEP_FILE                        Depth over time for each horizon'
write(0,*) '-time TIME_FILE                      Time file'
write(0,*) '-hf HF_FILE                          Surface heat flow over time'
write(0,*) '-tti TTI_FILE                        Time-Temperature Index of Maturity (TTI) over time'
write(0,*) '-closure FILE  T1[p] T2[p] T3[p]...  Closure temperature timing for each horizon'
write(0,*) '-temp:range FILE  T1-T2  T1-T2...    Closure temperature timing for each horizon'
write(0,*) '-format                              Print file formats and exit'
write(0,*)
if (printFileFormats) then
    write(0,*)
    write(0,*) 'File formats:'
    write(0,*) 'TEMP_FILE:'
    write(0,*) 'time_total(Ma)  timestep_interval(Ma)  number_timesteps'
    write(0,*) 'temp(t1,hor1) temp(t1,hor2) temp(t1,hor3)...'
    write(0,*) 'temp(t2,hor1) temp(t2,hor2) temp(t2,hor3)...'
    write(0,*) ':'
    write(0,*) ':'
    write(0,*)
    write(0,*) 'DEP_FILE:'
    write(0,*) 'time_total(Ma)  timestep_interval(Ma)  number_timesteps'
    write(0,*) 'dep(t1,hor1) dep(t1,hor2) dep(t1,hor3)...'
    write(0,*) 'dep(t2,hor1) dep(t2,hor2) dep(t2,hor3)...'
    write(0,*) ':'
    write(0,*) ':'
    write(0,*)
    write(0,*) 'TIME_FILE:'
    write(0,*) 'HF_FILE'
    write(0,*) 'CLOSURE_FILE:'
    write(0,*)
endif
stop
end subroutine
