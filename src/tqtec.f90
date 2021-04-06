module tqtec

! Inputs/outputs
character(len=512) :: input_file                  ! name of input file                              INFIL
character(len=8) :: input_mode                    ! how to read input parameters (user, file)
character(len=512) :: output_file                 ! name of output file                             OUTFIL

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
double precision, allocatable :: results(:,:,:)   ! temperature and depth for each timstep/depth    r1

end module


!==============================================================================!


program main

! C     CALCULATES THE ONE-DIMENSIONAL TRANSIENT THERMAL FIELD WITHIN
! C     AN AREA THAT UNDERGOES EPISODES OF:  BURIAL, EROSION, AND
! C     THRUSTING.  THIS PROGRAM CAN BE USED TO MONITOR POINTS THAT
! C     MOVE FROM THE UPPER TO LOWER (or v.v) PLATES.

! C     Q(V) = SURFACE HEAT FLOW AT EACH TIME STEP
! C     TTI(V,K) = TTI CALCULATED AT EACH TIME STEP


use tqtec

implicit none

integer :: ierr
integer :: i, np
double precision :: SURCON


! Initialize default model parameters
! Variable = value       ! Value in previous tqtec
nnodes = 1200            ! N=1200
nhorizons = 10           ! Hard-coded to 10
dt = 0.005d0             ! H1=0.005
dz = 0.05d0              ! K1=0.05
diffusivity = 32.0d0     ! D1=32.0


! Parse command line
! Matt's note: this is a totally new subroutine for tqtec (which I use in all my other programs)
! that allows better control over user input/output. Here, most of the model I/O is done via a
! control file, so gcmdln() is much simpler, only allowing specification of user or file inputs.
call gcmdln()
if (output_file.eq.'') then
    call usage('tqtec: output file must be defined')
endif


! Read control file or user input (formerly INPUT)
call read_model_parameters()


! Calculate model parameters and allocate arrays
nt_total = int(t_total/dt)
nt_output = int(t_output/dt)
temp_factor = diffusivity*dt/cond_base
r1 = diffusivity*dt/(dz*dz)
dtemp_wo_hp = (hf_surf-hp_surf*hp_dep)*dz/cond_base
allocate(depth_node(nhorizons))
allocate(conductivity(nnodes))
allocate(temp(nnodes))
allocate(hp(nnodes))
allocate(hf(nt_total))
allocate(results(nt_total,2,nhorizons))


! Set up action timing arrays (formerly: HIST)
call setup_action_arrays()


! Initialize the temperature, heat flow, and heat production at each node (formerly: INIT)
call initialize_thermal_parameters()


write(0,*) 'nnodes:        ',nnodes
write(0,*) 'nt_total:      ',nt_total
write(0,*) 'istep:         ',istep
write(0,*) 'dz:            ',dz
write(0,*) 'dt:            ',dt
write(0,*) 'r1:            ',r1
write(0,*) 't_total:       ',t_total
write(0,*) 't_output:      ',t_output
write(0,*) 'nt_output:     ',nt_output
write(0,*) 'nhorizons:     ',nhorizons
write(0,*) 'depth:         ',depth
write(0,*) 'depth_node:    ',depth_node
write(0,*) 'nlayers:       ',nlayers
write(0,*) 'layer(:,1):    ',layer(:,1)
write(0,*) 'layer(:,2):    ',layer(:,2)
write(0,*) 'layer(:,3):    ',layer(:,3)
write(0,*) 'diffusivity:   ',diffusivity
write(0,*) 'cond_base:     ',cond_base
write(0,*) 'temp_surf:     ',temp_surf
write(0,*) 'hp_surf:       ',hp_surf
write(0,*) 'hp_dep:        ',hp_dep
write(0,*) 'hf_surf:       ',hf_surf
write(0,*) 'hf_base:       ',hf_base
write(0,*) 'dtemp_wo_hp:   ',dtemp_wo_hp
write(0,*) 'dtemp_wo_hp:   ',dtemp_wo_hp
write(0,*) 'temp_factor:   ',temp_factor
write(0,*) 'temp_base_adj:  ',temp_base_adj
write(0,*) 'nburial:        ',nburial
write(0,*) 'burial_dat(:,1):',burial_dat(:,1)
write(0,*) 'burial_dat(:,2):',burial_dat(:,2)
write(0,*) 'burial_dat(:,3):',burial_dat(:,3)
write(0,*) 'burial_dat(:,4):',burial_dat(:,4)
write(0,*) 'nuplift:        ',nuplift
write(0,*) 'uplift_dat(:,1):',uplift_dat(:,1)
write(0,*) 'uplift_dat(:,2):',uplift_dat(:,2)
write(0,*) 'uplift_dat(:,3):',uplift_dat(:,3)
write(0,*) 'nthrust:        ',nthrust
write(0,*) 'thrust_dat(:,1):',thrust_dat(:,1)
write(0,*) 'thrust_dat(:,2):',thrust_dat(:,2)
write(0,*) 'thrust_dat(:,3):',thrust_dat(:,3)
write(0,*) 'thrust_dat(:,4):',thrust_dat(:,4)
write(0,*) 'thrust_dat(:,5):',thrust_dat(:,5)
write(0,*) 'nhfvars:        ',nhfvars
write(0,*) 'hfvar(:,1):     ',hfvar(:,1)
write(0,*) 'hfvar(:,2):     ',hfvar(:,2)


istep = 0
do while (istep.lt.nt_total)

    ! Update the adjusted temperature at the base of the model
! 5     W(3)=B(N)+BASGRAD(V+1)
    temp_base_adj = temp(nnodes) + bas_grad(istep+1)

    ! Calculate the updated temperatures at each node (the main finite difference procedure)
    ! (formerly: MAT and TRID)
    call update_temps(nnodes,ierr)
    if (ierr.ne.0) then
        write(0,*) 'tqtec: error in update_temps TRID algorithm'
        stop
    endif

    ! Increment the time step
    istep = istep + 1
    ! write(*,*) 'tqtec: working on cycle', istep

    ! Tectonic action!
    if (action(istep).eq.1) then
        call bury() ! (Formerly: BURIAL)
        ! print *,istep,'Burying...'
    elseif (action(istep).eq.2) then
        call erode() ! (Formerly: EROS)
        ! print *,istep,'Uplifting/eroding...'
    elseif (action(istep).ge.3) then
        if (int(thrust_dat(action(istep)-2,2)).eq.1) then
            call thrust_upperplate() ! (Formerly: THSTUP)
        elseif (int(thrust_dat(action(istep)-2,2)).eq.2) then
            call thrust_lowerplate() ! (Formerly: THSTLP)
        endif
    endif

!       Q(V)=(B(10)-B(5))/(5.0*H1)
    hf(istep) = (temp(10)-temp(5))/(5.0d0*dz)

!       SURCON=0.0
!       coninv = 0.0
!       DO 6, I=1,5
! C         coninv = coninv + (1.0/COND(I))
! 		 		 SURCON=SURCON+COND(I+4)
! 6     CONTINUE
!       SURCON=SURCON/5.0
    SURCON = 0.0d0
    do i = 1,5
        SURCON = SURCON + conductivity(i+4)
    enddo
    SURCON = SURCON/5.0d0
! C      SURCON = (1.0/coninv)/25.0

!       Q(V)=Q(V)*SURCON
    hf(istep) = hf(istep)*SURCON

!       DO 10 I=1,10
!          IARG=NINT(Y(I))
!          IF(IARG.EQ.0) THEN
!            R(V,1,I)=W(1)
!          ELSEIF(IARG.LT.0) THEN
!            R(V,1,I)=0.0
!          ELSEIF(IARG.GT.0) THEN
!            R(V,1,I)=B(IARG)
!          ENDIF
!          R(V,2,I)=Y(I)
! C         EMP = (R(V,1,I)-105)/10
! C         IF (V.EQ.1) THEN
! C           TTI(V,I)= II(2)*2**EMP
! C         ELSE
! C           TTI(V,I)= TTI(V-1,I) + II(2)*2**EMP
! C         ENDIF
! 10    CONTINUE
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

!       IF(V.GE.Q1)THEN
!         CALL OUTPUT(Q1,M1,OUTFIL)
!       ELSE
!         GOTO 5
!       ENDIF
    if (istep.ge.nt_total) then
        call output()
    endif
enddo



end


!--------------------------------------------------------------------------------------------------!


subroutine read_model_parameters()

use tqtec, only: input_mode, &
                 input_file, &
                 t_total, &
                 t_output, &
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
character(len=32) :: reply
character(len=512) :: input_line
logical :: inWhitespace


! write(0,*) 'read_model_parameters: starting'


if (input_mode.eq.'user') then

    write(*,*) 'Do you want to create a new data file? (Y/N)'
    read(*,*) reply
    if (reply.eq.'y'.or.reply.eq.'Y'.or.reply.eq.'1') then
        ! Create new data file with interactive input
        write(*,*) 'Name of input file to create?'
        read(*,*) input_file
        open(unit=9,file=input_file,status='unknown')
        write(9,*) trim(input_file)
    else
        write(*,*) 'Name of existing input file?'
        read(*,*) input_file
        call read_input_file()
        return
    endif

    write(*,*) 'Total time for model? (Ma)'
    read(*,*) t_total
    write(*,*) 'Time interval for output to be displayed? (Ma)'
    read(*,*) t_output
    t_output = 5.0d0 ! HARD-CODED

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

    write(9,*) t_total, t_output, temp_surf, hf_surf, cond_base, hp_surf, hp_dep


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

    write(9,*) nlayers
    if (nlayers.gt.0) then
        do i = 1,nlayers
            write(9,*) (layer(i,j),j=1,3)
        enddo
    endif


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
        write(*,*) 'Initial depth of point/horizon',i,'? (km)'
        read(*,*) depth(i)
    enddo

    write(9,*) (depth(i),i=1,nhorizons)


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

    write(9,*) nburial
    if (nburial.gt.0) then
        do i = 1,nburial
            write(9,*) (burial_dat(i,j),j=1,4)
        enddo
    endif


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

    write(9,*) nuplift
    if (nuplift.gt.0) then
        do i = 1,nuplift
            write(9,*) (uplift_dat(i,j),j=1,3)
        enddo
    endif


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

    write(9,*) nthrust
    if (nthrust.gt.0) then
        do i = 1,nthrust
            write(9,*) (thrust_dat(i,j),j=1,5)
        enddo
    endif


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

    write(9,*) nhfvars
    if (nhfvars.gt.0) then
        do i = 1,nhfvars
            write(9,*) (hfvar(i,j),j=1,2)
        enddo
    endif


    write(*,*) 'tqtec: input file "',trim(input_file),'" has been created'
    write(*,*) 'To re-use this file, run tqtec -f ',trim(input_file)


elseif (input_mode.eq.'file') then

    call read_input_file()

else
    write(0,*) 'tqtec: no input mode named "',trim(input_mode),'"'
    stop
endif

write(0,*) 'read_model_parameters: finished'


return
end subroutine


!--------------------------------------------------------------------------------------------------!

subroutine read_input_file()

use tqtec, only: input_mode, &
                 input_file, &
                 t_total, &
                 t_output, &
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
character(len=32) :: reply
character(len=512) :: input_line
logical :: inWhitespace

    ! Open the input file
    open(unit=8,file=input_file,iostat=ios)
    if (ios.ne.0) then
        write(0,*) 'read_model_parameters: something went wrong trying to open input file "', &
                   trim(input_file),'"'
        stop
    endif
    rewind(8)


    read(8,'(A)') input_line ! first line is name of file


    ! READ (8,110) Q1,M1,W(1),G1,C1,A1,B1
    read(8,'(A)') input_line ! second line contains first useful parameters
    read(input_line,*,iostat=ios) t_total, t_output, temp_surf, hf_surf, cond_base, hp_surf, &
                       hp_dep
    if (ios.ne.0) then
        read(input_line,*,iostat=ios) t_total, t_output, temp_surf, hf_surf, cond_base, hp_surf
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
    read(8,*) nburial
    if (allocated(burial_dat)) then
        deallocate(burial_dat)
    endif
    allocate(burial_dat(nburial,4))
    do i = 1,nburial
        read(8,*) (burial_dat(i,j),j=1,4)
    enddo


    ! Read uplift/erosion episodes
    read(8,*) nuplift
    if (allocated(uplift_dat)) then
        deallocate(uplift_dat)
    endif
    allocate(uplift_dat(nuplift,4))
    do i = 1,nuplift
        read(8,*) (uplift_dat(i,j),j=1,3)
    enddo


    ! Read thrust episodes
    read(8,*) nthrust
    if (allocated(thrust_dat)) then
        deallocate(thrust_dat)
    endif
    allocate(thrust_dat(nthrust,4))
    do i = 1,nthrust
        read(8,*) (thrust_dat(i,j),j=1,5)
    enddo


    ! ! Read basal heat flow variations
    ! read(8,*) nhfvars
    ! if (allocated(hfvar)) then
    !     deallocate(hfvar)
    ! endif
    ! allocate(hfvar(nhfvars,2))
    ! do i = 1,nhfvars
    !     read(8,*) (hfvar(i,j),j=1,2)
    ! enddo

return
end subroutine

!--------------------------------------------------------------------------------------------------!


subroutine setup_action_arrays()
!----
! Define arrays to control burial, erosion, and thrusting events
!----

use tqtec, only: nt_total, &
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


write(0,*) 'setup_action_arrays: starting'


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
    nstart = int(burial_dat(i,1)/dt)
    nduration = int(burial_dat(i,2)/dt)
    nthick = int(burial_dat(i,3)/dz)
    rate = dble(nthick)/dble(nduration)
    jbeg = nstart + 1
    jend = nstart + nduration
    ct = 0
    do j = jbeg,jend
        arg = dble(j-nstart)*rate - dble(ct)
        if (arg.ge.1.0d0) then
            action(j) = 1
            bcond(j) = burial_dat(i,4)
            ct = ct + 1
        endif
    enddo
enddo


! Uplift periods
do i = 1,nuplift
    nstart = int(uplift_dat(i,1)/dt)
    nduration = int(uplift_dat(i,2)/dt)
    nthick = int(uplift_dat(i,3)/dz)
    rate = dble(nthick)/dble(nduration)
    jbeg = nstart + 1
    jend = nstart + nduration
    ct = 0
    do j = jbeg,jend
        arg = dble(j-nstart)*rate - dble(ct)
        if (arg.ge.1.0d0) then
            action(j) = 2
            ct = ct + 1
        endif
    enddo
enddo


! Thrust events
do i = 1,nthrust
    nstart = int(thrust_dat(i,1)/dt)
    action(nstart) = i+2
    ! THTYPE(I): thrust_dat(i,2)
enddo


! Basal heat flow
! Initialize heat flow over time to be surface heat flow
hf_surf_var = hf_surf

! Timing of heat flow changes
do i = 1,nhfvars
    intqt(i) = int(hfvar(i,1)/dt)
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

! Propagate down to base
do j = 1,nt_total
    bas_grad(j) = (hf_surf_var(j)-hp_surf*hp_dep)*dz/cond_base
enddo


write(0,*) 'setup_action_arrays: finished'

return
end subroutine


!--------------------------------------------------------------------------------------------------!


subroutine initialize_thermal_parameters()
!----
! Calculate the steady state temperature at each node based on the surface heat flow, surface
! temperature, conductivity, and heat production
!----


use tqtec, only: nnodes, &
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
do i = 1,nnodes
    hp(i) = hp_surf*exp(-dble(i)*dz/hp_dep)
enddo


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


return
end subroutine


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!------------------------------- FINITE DIFFERENCE PROCEDURE --------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!


subroutine update_temps(nnodes,ierr)
!----
! A combination of old subroutines MAT (setting up the matrix equations for temperature) and TRID
! (solving the tridiagonal matrix equation for the new temperatures).
!----


use tqtec, only: r1, conductivity, temp, hp, temp_surf, temp_factor, temp_base_adj, dt, diffusivity

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



! write(0,*) 'update_temps: starting'


! Pre-calculate some terms for the thermal matrix
gam(1) = 1.0d0/conductivity(1)
bet(1) = 1.0d0/(gam(1)+gam(1))
do i = 2,nnodes
    gam(i) = 1.0d0/conductivity(i)
    bet(i) = 1.0d0/(gam(i)+gam(i-1))
enddo


! Load the thermal matrix (calculate local c, d, and e arrays that will be used later)
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


! Calculate temperature at each node
temp_new(1) = a(1,2)*temp(1) + a(1,3)*temp(2) + 2.0d0*a(1,1)*temp_surf + &
              diffusivity*dt*hp(1)/conductivity(1)
do i = 2,nnodes-1
    temp_new(i) = a(i,1)*temp(i-1) + a(i,2)*temp(i) + a(i,3)*temp(i+1) + &
                  diffusivity*dt*hp(i)/conductivity(i)
enddo
temp_new(nnodes) = a(nnodes,1)*temp(nnodes-1) + a(nnodes,2)*temp(nnodes) + &
               2.0d0*a(nnodes,3)*temp_base_adj + temp_factor*hp(nnodes)


! Update temperatures in global temperature array
temp = temp_new



! END SUBROUTINE MAT
! BEGIN SUBROUTINE TRID



! Initialize error flag
ierr = 0


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


! Backsubstitution
do l = 1,nnodes-2
    k = nnodes-2-l+1
    temp(k) = (temp(k)-d(k)*temp(k+1)-e(k)*temp(k+2))/c(k)
enddo


! write(0,*) 'update_temps: finished'


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

use tqtec, only: nnodes, &
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
! Erode and uplift the horizons by shifting physical parameters up node list and removing surface node
!----

use tqtec, only: nnodes, &
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
    if (depth_node(i).le.0) then
        depth_node(i) = 0
    endif
enddo

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine thrust_upperplate()
!----
! Generate a thrust fault, keeping the horizons in the upper plate
!----


! C     **************************************************************
! C     **************************************************************
! C     PROGRAM THSTUP.FOR:  THRUST FOR UPPER PLATE BOUNDARY
!       SUBROUTINE THSTUP(V)
! C     --------------------------------------------------------------

use tqtec, only: nnodes, &
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
integer :: i, k
integer :: thick_init, thrust_dep, thick_end, ierosion
double precision :: upl_conductivity(nnodes), upl_hp(nnodes), upl_temp(nnodes)


! Set the thrust number
k = action(istep) - 2

! Save thrust sheet parameters
thick_init = int(thrust_dat(k,3)/dz)  ! Thickness prior to thrusting, in nodes
thrust_dep = int(thrust_dat(k,4)/dz)  ! Depth of emplacement, in nodes
thick_end = int(thrust_dat(k,5)/dz)   ! Final thickness of thrust sheet, in nodes
if (thick_init.lt.thick_end) then
    write(0,*) 'thrust_upperplate: final thickness must be less than or equal to initial thickness'
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

! Update temperatures at tracked horizons
temp(1) = (temp_surf+temp(2))/2.0d0
temp(2) = (temp(1)+temp(2)+temp(3))/3.0d0
do i = 3,2*thick_init
    temp(i) = (temp(i-2)+temp(i-1)+temp(i)+temp(i+1)+temp(i+2))/5.0d0
enddo

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine thrust_lowerplate()
!----
! Generate a thrust fault, keeping the horizons in the upper plate
!----

use tqtec, only: nnodes, &
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
integer :: i, k
integer :: thick_init, thrust_dep, thick_end, ierosion, dnode
double precision :: upl_conductivity(nnodes), upl_hp(nnodes), upl_temp(nnodes)

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

! Update temperatures
temp(1) = (temp_surf+temp(2))/2.0d0
temp(2) = (temp(1)+temp(2)+temp(3))/3.0d0
do i = 3,2*thick_init
    temp(i) = (temp(i-2)+temp(i-1)+temp(i)+temp(i+1)+temp(i+2))/5.0d0
enddo

return
end subroutine


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
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
write(7,110) 0.0 ! II(8)
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
    write(7,130) depth(i)
enddo
130 format(F11.4)

return
end subroutine


!--------------------------------------------------------------------------------------------------!


subroutine gcmdln()
!----
! Parse tqtec command line arguments defining input/output modes and files
!----

use tqtec, only: input_mode, &
                 input_file, &
                 output_file

implicit none

! Local variables
character(len=512) arg
integer :: i, ios, narg


! Initialize control variables
ios = 0

! Initialize defaults
input_file = ''
input_mode = 'user'
output_file = ''


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

    elseif (arg.eq.'-interactive') then
        input_mode = 'user'

    elseif (arg.eq.'-o') then
        i = i + 1
        call get_command_argument(i,output_file,status=ios)

    else
        call usage('tqtec: no option '//trim(arg))
    endif

    9001 if (ios.ne.0) then
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
write(0,*) 'Usage: tqtec -interactive|-f INPUT_FILE  -o OUTPUT_FILE'
write(0,*)
write(0,*) '-interactive    User defines model parameters interactively'
write(0,*) '-f INPUT_FILE   Input model parameter file'
write(0,*) '-o OUTPUT_FILE  Output file'
write(0,*)
stop
end subroutine
