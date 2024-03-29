module readtqtec

! VARIABLE NAME -----------------------------------! DESCRIPTION ---------------------------------! PREVIOUS VARIABLE NAME

! Inputs/outputs
character(len=512) :: tqtec_output_file            ! tqtec output/readtqtec input file              INFILE
character(len=512) :: temp_file                    ! output temperature-timestep file               OUTFILE
character(len=512) :: dep_file                     ! output depth-timestep file                     OUTFILE
character(len=512) :: time_file                    ! output time-timestep file                      OUTFILE
character(len=512) :: hf_file                      ! output heat flow-timestep file                 OUTFILE
character(len=512) :: closure_file                 ! output closure temperature-depth file          OUTFILE
integer :: nclosure                                ! number of closure temperatures to track
double precision, allocatable :: closure_temps(:)  ! list of closure temperatures to track
logical, allocatable :: isPartialAnnealTemp(:)     ! flags for partial annealing
character(len=512) :: temp_range_file
integer :: ntemprange
double precision, allocatable :: temp_ranges(:,:)

! Finite difference parameters
integer :: nnodes                                 ! number of spatial nodes                         N
integer :: nt_total                               ! number of time steps                            Q1 (updated), II(5)
integer :: istep                                  ! current time step                               V
double precision :: dz                            ! node spacing (km)                               H1, II(1)
double precision :: dt                            ! time step interval (Ma)                         K1, II(2)
double precision :: r1                            ! finite difference time factor                   R1

! Timing
double precision :: t_total                       ! total model time (Ma)                           Q1 (initial), II(5)
double precision :: t_output                      ! time per output (Ma)                            M1 (initial)
integer :: nt_output                              ! time steps between outputs                      M1 (updated)

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
double precision, allocatable :: bas_grad(:)      !                                                 BASGRAD
integer, allocatable :: action(:)                 ! burial (1), erosion (2), or thrust (>=3)        P
double precision, allocatable :: bcond(:)         ! boundary condition magnitude                    BCOND

! Results array
double precision, allocatable :: results(:,:,:)   ! temperature and depth for each timestep/depth   r1

end module


!==================================================================================================!

program main
!----
! Read and print the output from tqtec
!----

use readtqtec

implicit none

! Local variables
integer :: i, j, k, l
character(len=1) :: dummy
character(len=32) :: fmt_string
logical :: ex, isClosed, isInRange
double precision :: xmin, closure_time, final_depth, temp1

! C     reads the output from program TQTec
!       CHARACTER INFILE*20, DUMMY*20, OUTFILE*20
!       REAL II
!       INTEGER Q1
!       DIMENSION S(7),E(7),C(7),A(7),DIFF(5000),PRODN(5000),
!      *   S2(9),E2(9),C2(9),A2(9),S3(13),E3(13),C3(13),A3(13)
!       COMMON /COM1/ II(10),Q(50000),R(50000,2,10),Y(10)


! Parse command line options
call gcmdln()


! Check whether tqtec output file exists
inquire(file=tqtec_output_file,exist=ex)
if (.not.ex) then
    call usage('readtqtec: no file found named "'//trim(tqtec_output_file)//'"')
endif


! Open the file and start reading
open(unit=8,file=tqtec_output_file,status='old')
!       OPEN (UNIT=8, FILE=INFILE)
!       REWIND 8

!       R1=1.99
! C     XMIN=0.0


! First line is just the file name
!       READ (8,'(A)') DUMMY
read(8,'(A)') dummy


! Read model parameters
!       DO 5 J=1,10
!          READ(8,*) II(J)
!          write(*,*) II(J)
! 5     CONTINUE
read(8,*) dz
read(8,*) dt
read(8,*) hp_surf
read(8,*) hp_dep
read(8,*) t_total
read(8,*) diffusivity
read(8,*) temp_surf ! temp_factor?
read(8,*) nhorizons
read(8,'(A)') dummy
read(8,'(A)') dummy


!       XMIN=II(5)
xmin = t_total

! Total number of timesteps
!       Q1=NINT(II(5)/(2*II(2)))
nt_total = int(t_total/(2.0d0*dt))
! 	  write (*,*) Q1


! Read heat flow at each timestep
!       DO 10 J=1,Q1
!          READ(8,*) Q(J)
! 10    CONTINUE
allocate(hf(nt_total))
do j = 1,nt_total
    read(8,*) hf(j)
enddo

! For each tracked horizon, read temperature and depth at each timestep
!       DO 25 K=1,10
!          DO 26 J=1,2
!             DO 27 I=1,Q1
!                READ(8,*) R(I,J,K)
! 27          CONTINUE
! 26       CONTINUE
! 25    CONTINUE
allocate(results(nt_total,2,nhorizons))
do k = 1,nhorizons
    do j = 1,2
        do i = 1,nt_total
            read(8,*) results(i,j,k)
        enddo
    enddo
enddo

! C      DO 30 K=1,5
! C         DO 30 I=1,Q1
! C            READ(8,120) TTI(I,K)
! C 30    CONTINUE

! Read final depth of each horizon
!       DO 40 I=1,10
!          READ(8,*) Y(I)
! 40    CONTINUE
allocate(depth_node(nhorizons))
do i = 1,nhorizons
    read(8,*) depth_node(i)
enddo


! 100   FORMAT(A20)
! 110   FORMAT(F7.3)
! 115   FORMAT(F6.2)
! 120   FORMAT(F7.1)
! 130   FORMAT(F11.4)
! 170   FORMAT(I1)

! Close the TQTec output file
CLOSE(8)


! C
! C     CHOOSE PLOT DATA DESIRED
! C
! 500   WRITE(*,*)'WHICH PLOT DATA DO YOU WANT? (ONE DIGIT ONLY)'
!       WRITE(*,*)'TEMPERATURE =(1)'
!       WRITE(*,*)'DEPTH       =(2)'
! C      WRITE(*,*)'TTI         =(3)'
! C      WRITE(*,*)'PROD. TYPE 1=(4)'
! C      WRITE(*,*)'PROD. TYPE 2=(5)'
! C      WRITE(*,*)'PROD. TYPE 3=(6)'
!       WRITE(*,*)'SURFACE Q   =(3)'
! 	  WRITE(*,*) 'TIME@OUTPUT  = (4)'
!       WRITE(*,*)'QUIT?       =(9)'
!       READ(*,170) IPLT
!       IF (IPLT.EQ.1) GOTO 510
!       IF (IPLT.EQ.2) GOTO 520
!       IF (IPLT.EQ.3) GOTO 570
! 	  IF (IPLT.EQ.4) GOTO 580
! C      IF (IPLT.EQ.4) GOTO 540
! C      IF (IPLT.EQ.5) GOTO 550
! C      IF (IPLT.EQ.6) GOTO 560
! C      IF (IPLT.EQ.7) GOTO 570
!       IF (IPLT.EQ.9) GOTO 590


! C
! C     WRITE TEMPERATURE PLOT DATA
! C
! 510   WRITE(*,*) 'TYPE NAME OF OUTPUT FILE FOR TEMP Data'
!       READ(*,100) OUTFILE
!       OPEN (UNIT=9,FILE=OUTFILE)
!       WRITE(9,175) -XMIN,II(2)*2,Q1
!       DO 50 L=1,Q1
! 		WRITE(9,160)(R(L,1,I),I=1,10)
! 50    CONTINUE
!       CLOSE(9)
!       GOTO 500
if (temp_file.ne.'') then
    open(unit=9,file=temp_file,status='unknown')
    write(9,'(2F10.3,I6)') -xmin,2.0d0*dt,nt_total          ! xmin=time before present (Ma)
    write(fmt_string,'("("I6,"F8.3)")') nhorizons           ! 2*dt=sample rate (Ma) (because saved every 2 points from tqtec run)
    do l = 1,nt_total
        write(9,fmt_string) (results(l,1,i),i=1,nhorizons)
    enddo
    close(9)
endif


! C
! C     WRITE DEPTH PLOT DATA
! C
! 520   WRITE(*,*)'TYPE NAME OF OUTPUT FILE FOR DEPTH Data'
!       READ(*,100) OUTFILE
!       OPEN (UNIT=9, FILE=OUTFILE)
!       WRITE(9,175) -XMIN,II(2)*2,Q1
!       DO 60 L=1,Q1
!           WRITE(9,160)(-1*R(L,2,I)*II(1),I=1,10)
! 60    CONTINUE
!       CLOSE (9)
!       GOTO 500
if (dep_file.ne.'') then
    open(unit=9,file=dep_file,status='unknown')
    write(9,'(2F10.3,I6)') -xmin,2.0d0*dt,nt_total 
    write(fmt_string,'("("I6,"F8.3)")') nhorizons
    do l = 1,nt_total
        write(9,fmt_string) (-1.0d0*results(l,2,i)*dz,i=1,nhorizons)
    enddo
    close(9)
endif

! C
! C     WRITE TTI PLOT DATA
! C
! C530   WRITE(*,*)'TYPE NAME OF OUTPUT FILE FOR TTI PLOTS'
! C      READ(*,100) OUTFILE
! C      OPEN (UNIT=9,FILE=OUTFILE)
! C         WRITE(9,175) -XMIN,II(2)*2,Q1,LOG10(1.),LOG10(100000.)
! C      DO 70 L=1,5
! C         DO 71 M=1,Q1
! C            IF (TTI(M,L).LE.1.0) THEN
! C               TTI(M,L)=0.0
! C            ELSE
! C               TTI(M,L)=LOG10(TTI(M,L))
! C            ENDIF
! C71       CONTINUE
! C         DO 75 K=10,Q1,10
! C75          WRITE(9,165)(TTI(I,L),I=K-9,K)
! C         IF (K.NE.Q1.AND.K.NE.(Q1+10)) THEN
! C            KNEW=ABS(Q1-K)
! C            WRITE(9,165)(TTI(I,L),I=Q1-KNEW+1,Q1)
! C         ENDIF
! C70    CONTINUE
! C      CLOSE(9)
! C      GOTO 500
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

! C
! C     WRITE SURFACE HEAT FLOW DATA
! C
! 570   WRITE(*,*) 'TYPE NAME OF OUTPUT FILE FOR SURF. Q Data'
!       READ(*,100) OUTFILE
!       OPEN (UNIT=9,FILE=OUTFILE)
!       WRITE(9,175) -XMIN,II(2)*2,Q1
! 	  DO 575 I=1,Q1
! 575      WRITE(9,160) Q(I)
! 	  CLOSE(9)
!       GOTO 500
! C      Q(2)=(Q(1)+Q(2)+Q(3))/3.0
! C      Q(Q1-1)=(Q(Q1)+Q(Q1-1)+Q(Q1-2))/3.0
! C      DO 571 I=3,Q1-2
! C         Q(I)=(Q(I-2)+Q(I-1)+Q(I)+Q(I+1)+Q(I+2))/5
! C 571   CONTINUE
if (hf_file.ne.'') then
    open(unit=9,file=hf_file,status='unknown')
    write(9,'(2F10.3,I6)') -xmin,2.0d0*dt,nt_total
    do i = 1,nt_total
        write(9,*) hf(i)
    enddo
    close(9)
endif

! C
! C     WRITE TIME DATA
! C
! 580   WRITE(*,*)'TYPE NAME OF OUTPUT FILE FOR TIME Data'
!       READ(*,100) OUTFILE
!       OPEN (UNIT=9, FILE=OUTFILE)
!       WRITE(9,175) -XMIN,II(2)*2,Q1
!       DO 75 L=1,Q1
!           WRITE(9,160)(L*II(2)*2,I=1,10)
! 75    CONTINUE
!       CLOSE (9)
!       GOTO 500
if (time_file.ne.'') then
    open(unit=9,file=time_file,status='unknown')
    write(9,'(2F10.3,I6)') -xmin,2.0d0*dt,nt_total
    write(fmt_string,'("("I6,"F8.3)")') nhorizons
    do l = 1,nt_total
        write(9,fmt_string) (l*dt*2.0d0,i=1,nhorizons)
    enddo
    close(9)
endif


!----
! Write closure temperature timing
!----
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


! C
! 150   FORMAT(I2,2F7.1,I7)
! 160   FORMAT(10F8.3)
! 165   FORMAT(10F7.3)
! 175   FORMAT(2F10.3,I6)
! 590   CONTINUE
!       STOP
!       END
end

!--------------------------------------------------------------------------------------------------!


subroutine gcmdln()

use readtqtec, only: tqtec_output_file, &
                     temp_file, &
                     dep_file, &
                     time_file, &
                     hf_file, &
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
