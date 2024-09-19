!--------------------------------------------------------------------------------------------------!
! MinAge                                                                                           !
!                                                                                                  !
! Authors:                                                                                         !
!     - Matt Herman (current Modern Fortran version, i.e., what you are looking at right now!)     !
!     - Matt Legg, Rachel Piotraschke (PSU MS theses)                                              !
!                                                                                                  !
! MinAge (Mineral Age) calculates geochronological observables for a temperature-(depth-)time      !
! history. Currently, minage supports calculating:                                                 !
!     - Apatite Fission Track Ages/Distributions                                                   !
!     - Apatite (Uranium-Thorium)/Helium                                                           !
!                                                                                                  !
! minage is designed to work with TQTec and readTQTec output files, but will work with any         !
! temperature history in the format:                                                               !
!                                                                                                  !
!     total_time(Ma) timestep(Ma) nsteps                                                           !
!     temp_1                                                                                       !
!     temp_2                                                                                       !
!     temp_3                                                                                       !
!       :                                                                                          !
!     temp_nsteps                                                                                  !
!--------------------------------------------------------------------------------------------------!


module minage


    ! Input files
    character(len=512) :: readtqtec_temp_file
    character(len=512) :: readtqtec_dep_file

    ! Output files
    character(len=512) :: ahe_file
    logical :: isOutputDefined

    ! Number of horizons to track
    integer :: nhorizons
    integer, parameter :: nhorizons_max = 100

    ! Temperature-depth time series variables
    double precision :: time_ma
    double precision :: dt_ma
    integer :: ntimes
    double precision, allocatable :: temp_celsius_array(:,:)
    double precision, allocatable :: dep_km_array(:,:)

    ! Apatite fission track variables
    character(len=512) :: aft_file
    character(len=512) :: aft_history_file
    double precision :: aft_history_dt
    integer :: aft_n0
    double precision, allocatable :: aft_len0(:)
    logical :: doSegmentationCarlson1990
    logical :: doEtchingUserBiasCorrectionWillett1997

    ! Apatite (U-Th)/He variables
    integer :: ahe_nradius
    double precision, allocatable :: ahe_radius(:)
    double precision :: ahe_beta
    double precision :: ahe_dt_var
    integer :: ahe_nnodes
    double precision :: ahe_taumax
    integer :: ahe_verbosity



end module


!==================================================================================================!


program main


use minage, only:        &
    readtqtec_temp_file, &
    readtqtec_dep_file,  &
    aft_file,            &
    aft_history_file,    &
    ahe_file,            &
    isOutputDefined


implicit none




! Parse command line arguments
call gcmdln()




! Input/output error checks

! Is input temperature history defined?
if (readtqtec_temp_file.eq.'') then
    call usage('minage: no readtqtec output temperature file defined')
endif

! Is input depth history defined?
if (readtqtec_dep_file.eq.'') then
    write(0,*) 'minage: no readtqtec output depth file defined - setting all depths to 0 km'
endif

! Is an output defined?
if (.not.isOutputDefined) then
    call usage('minage: no output defined')
endif




! Read temperature history
call read_temp_history()



! Read depth history
call read_dep_history()





!************************************!
!********** CALCULATE AGES **********!
!************************************!

! Apatite fission track (variables and descriptions in fission_track_module.f90)
if (aft_file.ne.''.or.aft_history_file.ne.'') then
    call calc_aft_ages()
endif


! Apatite (U-Th)/He (variables and descriptions in radiogenic_helium_module.f90)
if (ahe_file.ne.'') then
    call calc_ahe_ages()
endif



end program





!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!---------------------------------- INPUT SUBROUTINES ---------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!



subroutine read_temp_history()
!----
! Program readtqtec reads output from TQTec and prints temperature [C] over time [Ma] for each
! user-specified horizon (-temp option). Read this temperature-time history for each horizon into
! an array.
!----

use minage, only: readtqtec_temp_file, &
                  nhorizons, &
                  nhorizons_max, &
                  time_ma, &
                  dt_ma, &
                  ntimes, &
                  temp_celsius_array

implicit none

! Local variables
integer :: i
integer :: ihorizon
integer :: ios
character(len=512) :: line
real :: dum(nhorizons_max)


! Open file for reading
open(unit=21,file=readtqtec_temp_file,status='old')


! First line: total time [Ma], sample rate [Ma], and number of temperature-time points
read(21,*,iostat=ios) time_ma, dt_ma, ntimes
if (ios.ne.0) then
    write(0,*) 'minage: error reading first line of tqtec output temperature file'
    call error_exit(1)
endif


! Extract number of horizons tracked in the thermal model from first line of temperatures
nhorizons = 1
read(21,'(A)',iostat=ios) line
if (ios.ne.0) then
    write(0,*) 'minage: error reading second line of tqtec output temperature file'
    call error_exit(1)
endif
do i = 1,nhorizons_max
    read(line,*,iostat=ios) dum(1:i)
    if (ios.ne.0) then
        nhorizons = i-1
        exit
    endif
enddo


! Step back one line in temp file to read all temperatures
backspace(21)


! Allocate and initialize the temperature array
allocate(temp_celsius_array(nhorizons,ntimes))
temp_celsius_array = 0.0d0


! Fill the temperature array
do i = 1,ntimes
    read(21,*,iostat=ios) (temp_celsius_array(ihorizon,i),ihorizon=1,nhorizons)
    if (ios.ne.0) then
        write(0,*) 'minage: error reading line',i,' of readtqtec output temperature file'
        call error_exit(1)
    endif
enddo


! Close the temperature file
close(21)


return
end subroutine



!--------------------------------------------------------------------------------------------------!



subroutine read_dep_history()
!----
! Program readtqtec reads output from TQTec and prints depth [C] over time [Ma] for each
! user-specified horizon (-dep option). Read this depth-time history for each horizon into an
! array.
!----

use minage, only: readtqtec_dep_file, &
                  nhorizons, &
                  nhorizons_max, &
                  time_ma, &
                  dt_ma, &
                  ntimes, &
                  dep_km_array

implicit none

! Local variables
integer :: i
integer :: ihorizon
integer :: nhorizons_dep
integer :: ntimes_dep
integer :: ios
character(len=512) :: line
real :: dum(nhorizons_max)



! If no depth file provided, keep things simple
if (readtqtec_dep_file.eq.'') then
    allocate(dep_km_array(nhorizons,ntimes))
    dep_km_array = 0.0d0
    return
endif


! Otherwise...

! Open file for reading
open(unit=20,file=readtqtec_dep_file,status='old')


! First line: total time [Ma], sample rate [Ma], and number of depth-time points
read(20,*,iostat=ios) time_ma, dt_ma, ntimes_dep
if (ios.ne.0) then
    write(0,*) 'minage: error reading first line of tqtec output depth file'
    call error_exit(1)
endif
if (ntimes_dep.ne.ntimes) then
    call usage('minage: number of timesteps differs between temp and dep files')
endif


! Extract number of horizons tracked in the thermal model from first line of temperatures
nhorizons_dep = 1
read(20,'(A)',iostat=ios) line
if (ios.ne.0) then
    write(0,*) 'minage: error reading second line of tqtec output depth file'
    call error_exit(1)
endif
do i = 1,nhorizons_max
    read(line,*,iostat=ios) dum(1:i)
    if (ios.ne.0) then
        nhorizons_dep = i-1
        exit
    endif
enddo
if (nhorizons_dep.ne.nhorizons) then
    call usage('minage: number of tracked horizons differs between temp and dep files')
endif


! Step back one line in depth file to read all depths
backspace(20)


! Allocate and initialize the depth array
allocate(dep_km_array(nhorizons,ntimes))
dep_km_array = 0.0d0


! Fill the depth array
do i = 1,ntimes
    read(20,*,iostat=ios) (dep_km_array(ihorizon,i),ihorizon=1,nhorizons)
    if (ios.ne.0) then
        write(0,*) 'minage: error reading line',i,' of readtqtec output depth file'
        call error_exit(1)
    endif
enddo


! Close the depth file
close(20)


return
end subroutine






!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!----------------------------- AGE CALCULATION SUBROUTINES ----------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!




subroutine calc_aft_ages()
!----
! Calculate apatite fission track ages from temperature-depth history. Subroutines and variables
! for fission track calculations can be found in fission_track_module.f90
!----


use minage, only:                          &
    aft_file,                              &
    aft_history_file,                      &
    aft_history_dt,                        &
    nhorizons,                             &
    dt_ma,                                 &
    ntimes,                                &
    temp_celsius_array,                    &
    dep_km_array,                          &
    aft_n0,                                &
    aft_len0,                              &
    doSegmentationCarlson1990,             &
    doEtchingUserBiasCorrectionWillett1997


use fission_track, only:                       &
    generate_fts_carlson_1990,                 &
    segment_fts_carlson_1990,                  &
    correct_fts_etching_userbias_willett_1997, &
    calc_ft_retention_age,                     &
    calc_ft_age,                               &
    load_histogram


implicit none


! Local variables
integer :: i
integer :: idt
integer :: itm
integer :: ibin
integer :: nbins
integer :: ihorizon
integer, allocatable :: aft_hist(:,:)
integer, allocatable :: aft_hist_corr(:,:)
double precision :: aft_len0_mean
double precision :: len
double precision :: xmin
double precision :: binwid
double precision, allocatable :: aft_len(:,:)
double precision, allocatable :: aft_age(:)
double precision, allocatable :: aft_retention_age(:)
character(len=32) :: fmt_string
logical :: isOpen



write(*,*) 'minage: calculating apatite fission track ages'


!----
!   Initialize arrays
!----
! Lengths of fission tracks generated during a fission event
if (aft_n0.eq.0) then    ! if not set on the command line, use default values (Green et al., 1986)
    aft_n0 = 20
    allocate(aft_len0(aft_n0))
    aft_len0(1:5)   = 15.0d0
    aft_len0(6:12)  = 16.0d0
    aft_len0(13:19) = 17.0d0
    aft_len0(20)    = 18.0d0
endif

! Calculate average initial fission track length
aft_len0_mean = 0.0d0
do i = 1,aft_n0
    aft_len0_mean = aft_len0_mean + aft_len0(i)
enddo
aft_len0_mean = aft_len0_mean/dble(aft_n0)


! Fission track length array
allocate(aft_len(aft_n0,ntimes))


! Histogram parameters (THESE COULD BE COMMAND LINE VARIABLES...)
xmin = 1.0d0
binwid = 1.0d0
nbins = 20
allocate(aft_hist(nbins,nhorizons))
allocate(aft_hist_corr(nbins,nhorizons))

! Fission track age arrays
allocate(aft_retention_age(nhorizons))
allocate(aft_age(nhorizons))



!*****************************************************************************!
! Calculate fission track distribution and ages at end of temperature history !
!*****************************************************************************!
if (aft_file.ne.'') then

    open(unit=22,file=aft_file,status='unknown')

    do ihorizon = 1,nhorizons ! for each horizon...

        ! Calculate fission track lengths for the temperature history -> aft_len
        call generate_fts_carlson_1990(     &
            aft_n0,                         &
            aft_len0,                       &
            ntimes,                         &
            temp_celsius_array(ihorizon,:), &
            dep_km_array(ihorizon,:),       &
            dt_ma,                          &
            'green-et-al-1986',             &
            aft_len                         &
        )


        ! Load the fission track lengths into a histogram -> aft_hist
        call load_histogram(     &
            aft_n0*ntimes,       &
            aft_len,             &
            binwid,              &
            xmin,                &
            nbins,               &
            aft_hist(:,ihorizon) &
        )
        aft_hist_corr(:,ihorizon) = aft_hist(:,ihorizon)


        ! Correct for fission track segmentation using Carlson (1990) method
        if (doSegmentationCarlson1990) then
            call segment_fts_carlson_1990( &
                nbins,                     &
                binwid,                    &
                xmin,                      &
                aft_hist_corr(:,ihorizon)  &
            )
        endif


        ! Correct for user bias and etching using Willett (1997) method
        if (doEtchingUserBiasCorrectionWillett1997) then
            call correct_fts_etching_userbias_willett_1997( &
                nbins,                                      &
                binwid,                                     &
                xmin,                                       &
                aft_len0_mean,                              &
                aft_hist_corr(:,ihorizon)                   &
            )
        endif


        ! Calculate fission track retention age (based on number of tracks)
        call calc_ft_retention_age(     &
            nbins,                      &
            aft_hist_corr(:,ihorizon),  &
            aft_n0,                     &
            dt_ma,                      &
            aft_retention_age(ihorizon) &
        )


        ! Calculate fission track age (based on number of tracks and average track length)
        call calc_ft_age(              &
            nbins,                     &
            binwid,                    &
            xmin,                      &
            aft_hist_corr(:,ihorizon), &
            aft_n0,                    &
            aft_len0_mean,             &
            dt_ma,                     &
            aft_age(ihorizon)          &
        )


        ! ! Save temperature at time of track-count-based (retention) age
        ! if (.not.allocated(temp_ft_retention_age_corr)) then
        !     allocate(temp_ft_retention_age_corr(nhorizons))
        ! endif
        ! if (time_ma.lt.0.0d0) then
        !     itemp = nint(-time_ma/dt_ma - ft_retention_age_corr/dt_ma)
        ! else
        !     itemp = nint(time_ma/dt_ma - ft_retention_age_corr/dt_ma)
        ! endif
        ! if (itemp.lt.1) then
        !     itemp = 1
        ! elseif (itemp.gt.ntimes) then
        !     itemp = ntimes
        ! endif
        ! temp_ft_retention_age_corr(ihorizon) = temp_celsius_array(ihorizon,itemp)


    enddo   ! end of loop over horizons


    ! File header
    write(22,'(A)')   '# Apatite fission track ages (Ma)'
    write(22,'(A)')   '#'
    write(22,'(A)')   '# CORRECTIONS:'
    write(22,'(A,L)') '# doSegmentationCarlson1990=',doSegmentationCarlson1990
    write(22,'(A,L)') '# doEtchingUserBiasCorrectionWillett1997=',doEtchingUserBiasCorrectionWillett1997
    write(22,'(A)')   '#'

    ! Retention ages
    write(22,'(A)')   '# RETENTION AGE'
    write(fmt_string,'("(6X,"I6,"F8.3)")') nhorizons
    write(22,fmt_string) (aft_retention_age(ihorizon),ihorizon=1,nhorizons)
    write(22,'(A)')   '#'

    ! Fission track ages
    write(22,'(A)')   '# FISSION TRACK AGE'
    write(fmt_string,'("(6X,"I6,"F8.3)")') nhorizons
    write(22,fmt_string) (aft_age(ihorizon),ihorizon=1,nhorizons)
    write(22,'(A)')   '#'

    ! write(22,*) 'Temperatures (C) at time of track-count-based (retention) ages'
    ! write(fmt_string,'("(6X,"I6,"F8.2)")') nhorizons
    ! write(22,fmt_string) (temp_ft_retention_age_corr(ihorizon),ihorizon=1,nhorizons)
    ! write(22,*)

    ! Track length histograms
    write(22,'(A)') '# TRACK LENGTH HISTOGRAMS'
    write(22,'(A)') '# Final (corrected) track lengths'
    write(fmt_string,'("(F6.1,"I6,"I8)")') nhorizons
    do ibin = 1,nbins
        len = xmin+dble(ibin-1)*binwid
        write(22,fmt_string) len, (aft_hist_corr(ibin,ihorizon),ihorizon=1,nhorizons)
    enddo
    write(22,'(A)')   '#'
    write(22,'(A)') '# Initial (uncorrected) track lengths'
    do ibin = 1,nbins
        len = xmin+dble(ibin-1)*binwid
        write(22,fmt_string) len, (aft_hist(ibin,ihorizon),ihorizon=1,nhorizons)
    enddo

    ! Close file
    close(22)
endif




!*************************************************************************!
! Calculate fission track distributions and ages over temperature history !
!*************************************************************************!
if (aft_history_file.ne.'') then

    open(unit=23,file=aft_history_file,status='unknown')


    idt = int(aft_history_dt/dt_ma) ! Interval in timesteps

    do itm = idt,ntimes,idt

        do ihorizon = 1,nhorizons

            ! Calculate fission track lengths at the time of interest -> aft_len
            call generate_fts_carlson_1990(         &
                aft_n0,                             &
                aft_len0,                           &
                itm,                                &
                temp_celsius_array(ihorizon,1:itm), &
                dep_km_array(ihorizon,1:itm),       &
                dt_ma,                              &
                'green-et-al-1986',                 &
                aft_len(1:aft_n0,1:itm)             &
            )

            ! Load the fission track lengths into a histogram -> aft_history_hist
            call load_histogram(         &
                aft_n0*itm,              &
                aft_len(1:aft_n0,1:itm), &
                binwid,                  &
                xmin,                    &
                nbins,                   &
                aft_hist(:,ihorizon)     &
            )

            ! Calculate fission track retention age (based on number of tracks)
            call calc_ft_retention_age(     &
                nbins,                      &
                aft_hist(:,ihorizon),       &
                aft_n0,                     &
                dt_ma,                      &
                aft_retention_age(ihorizon) &
            )

            ! Calculate fission track age (based on number of tracks and average track length)
            call calc_ft_age(              &
                nbins,                     &
                binwid,                    &
                xmin,                      &
                aft_hist(:,ihorizon),      &
                aft_n0,                    &
                aft_len0_mean,             &
                dt_ma,                     &
                aft_age(ihorizon)          &
            )

        enddo

        ! Write results to file
        write(23,'(A,F8.3)') '# MODEL_TIME',dble(itm)*dt_ma
        write(fmt_string,'("(A,2X,"I6,"F8.3)")') nhorizons
        write(23,fmt_string) '# RETENTION_AGE',(aft_retention_age(ihorizon),ihorizon=1,nhorizons)
        write(fmt_string,'("(A,2X,"I6,"F8.3)")') nhorizons
        write(23,fmt_string) '# AFT_AGE      ',(aft_age(ihorizon),ihorizon=1,nhorizons)
        write(fmt_string,'("(4X,F6.1,7X,"I6,"I8)")') nhorizons
        do ibin = 1,nbins
            len = xmin+dble(ibin-1)*binwid
            write(23,fmt_string) len, (aft_hist(ibin,ihorizon),ihorizon=1,nhorizons)
        enddo
    enddo

    ! Close file
    inquire(unit=23,opened=isOpen)
    if (isOpen) then
        close(23)
    endif

endif



! Free array memory
if (allocated(aft_len)) then
    deallocate(aft_len)
endif
if (allocated(aft_hist)) then
    deallocate(aft_hist)
endif
if (allocated(aft_hist_corr)) then
    deallocate(aft_hist_corr)
endif
if (allocated(aft_age)) then
    deallocate(aft_age)
endif
if (allocated(aft_retention_age)) then
    deallocate(aft_retention_age)
endif


return

end subroutine




!--------------------------------------------------------------------------------------------------!




subroutine calc_ahe_ages()
!----
! Calculate apatite (U-Th)/He ages from a temperature-depth history. Subroutines and variables
! for apatite-helium calculations can be found in radiogenic_helium_module.f90
!----


use minage, only: &
    ahe_file, &
    nhorizons, &
    dt_ma, &
    ntimes, &
    temp_celsius_array, &
    dep_km_array, &
    ahe_nradius, &
    ahe_radius, &
    ahe_nnodes, &
    ahe_beta, &
    ahe_dt_var, &
    ahe_taumax, &
    ahe_verbosity

use radiogenic_helium, only: &
    calc_apatite_he_age


implicit none


! Local variables
integer :: ihorizon
integer :: iradius
double precision :: radius_microns
double precision, allocatable :: ahe_age(:,:)
character(len=32) :: fmt_string



write(*,*) 'minage: calculating apatite (U-Th)/He ages'


! Open apatite-helium output file
open(unit=23,file=ahe_file,status='unknown')



! Apatite grain radii
if (ahe_nradius.eq.0) then    ! if not set on the command line, use default values
    ahe_nradius = 7
    allocate(ahe_radius(ahe_nradius))
    ahe_radius(1) = 20.0d0
    ahe_radius(2) = 40.0d0
    ahe_radius(3) = 60.0d0
    ahe_radius(4) = 80.0d0
    ahe_radius(5) = 100.0d0
    ahe_radius(6) = 120.0d0
    ahe_radius(7) = 140.0d0
endif



! Array for grain-size-based ages in each tracked horizon
if (.not.allocated(ahe_age)) then
    allocate(ahe_age(nhorizons,ahe_nradius))
endif
ahe_age = 0.0d0


! For each horizon...
do ihorizon = 1,nhorizons


    ! For each grain radius...
    do iradius = 1,ahe_nradius


        if (ahe_verbosity.ge.1) then
            write(*,2301) 'minage: working on horizon',ihorizon,'of',nhorizons, &
                        ': grain size',iradius,'of',ahe_nradius
            2301 format(X,2(A,I6,X,A,I6))
        endif


        ! Calculate (U-Th)/He age for this temperature-depth-time history and grain radius
        radius_microns = ahe_radius(iradius)
        call calc_apatite_he_age(           &
            ntimes,                         &
            temp_celsius_array(ihorizon,:), &
            dep_km_array(ihorizon,:),       &
            dt_ma,                          &
            ahe_dt_var,                     &
            radius_microns,                 &
            ahe_nnodes,                     &
            ahe_beta,                       &
            ahe_taumax,                     &
            ahe_verbosity,                  &
            ahe_age(ihorizon,iradius)       &
        )


    enddo

enddo


! Print results to apatite helium file
write(23,*) '# Apatite (U-Th)/He ages (Ma)'
write(23,*) '# Radius(microns) Age(Ma)'
write(fmt_string,'("(F12.1,"I6,"F12.3)")') nhorizons
do iradius = 1,ahe_nradius
    write(23,fmt_string) ahe_radius(iradius), &
                         (ahe_age(ihorizon,iradius),ihorizon=1,nhorizons)
enddo


! Close file
close(23)


return

end subroutine



!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!------------------------------------ INPUT/OUTPUT SUBROUTINES ------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!



subroutine gcmdln()
!----
! Parse command line
!----

use minage, only: &
    readtqtec_temp_file, &
    readtqtec_dep_file, &
    ahe_file, &
    isOutputDefined, &
    aft_file, &
    aft_history_file, &
    aft_history_dt, &
    aft_n0, &
    aft_len0, &
    doSegmentationCarlson1990, &
    doEtchingUserBiasCorrectionWillett1997, &
    ahe_nradius, &
    ahe_radius, &
    ahe_beta, &
    ahe_dt_var, &
    ahe_nnodes, &
    ahe_taumax, &
    ahe_verbosity


implicit none

integer :: i
integer :: j
integer :: narg
integer :: ios
character(len=512) :: tag
logical :: fileExists


! Initialize variables
readtqtec_temp_file = ''
readtqtec_dep_file = ''
ahe_file = ''
isOutputDefined = .false.


! Apatite fission track variables
aft_file = ''
aft_history_file = ''
aft_history_dt = 1.0d0
aft_n0 = 0
doSegmentationCarlson1990 = .true.
doEtchingUserBiasCorrectionWillett1997 = .true.

! Apatite (U-Th)/He variables
ahe_nradius = 0
ahe_dt_var = 0.1d0
ahe_nnodes = 102
ahe_beta = 0.85d0
ahe_taumax = 0.40d0
ahe_verbosity = 0


! Count number of command line arguments
narg = command_argument_count()
if (narg.eq.0) then
    call usage('')
endif


! Read and parse command line
i = 1
do while (i.le.narg)

    call get_command_argument(i,tag)

    ! Input files
    if (trim(tag).eq.'-temp') then
        i = i + 1
        call get_command_argument(i,readtqtec_temp_file,status=ios)

    elseif (trim(tag).eq.'-dep') then
        i = i + 1
        call get_command_argument(i,readtqtec_dep_file,status=ios)
        inquire(file=readtqtec_dep_file,exist=fileExists)
        if (.not.fileExists) then
            write(0,*) 'minage: could not find depth file ',trim(readtqtec_dep_file)
            readtqtec_dep_file = ''
        endif

    ! Output files
    elseif (trim(tag).eq.'-ahe') then
        i = i + 1
        call get_command_argument(i,ahe_file,status=ios)
        isOutputDefined = .true.


    !*************************************************!
    ! Apatite fission track age calculation variables !
    !*************************************************!
    elseif (trim(tag).eq.'-aft') then
        i = i + 1
        call get_command_argument(i,aft_file,status=ios)
        isOutputDefined = .true.

    elseif (trim(tag).eq.'-aft:history') then
        i = i + 1
        call get_command_argument(i,tag,status=ios)
        if (ios.ne.0) then
            write(0,*) 'minage: reached end of line before -aft:history arguments read'
            call error_exit(1)
        endif
        read(tag,*,iostat=ios) aft_history_dt
        if (ios.ne.0) then
            write(0,*) 'minage: tried to read -aft:history DT but could not parse "',trim(tag),'"'
            call error_exit(1)
        endif
        i = i + 1
        call get_command_argument(i,aft_history_file)
        isOutputDefined = .true.

    elseif (trim(tag).eq.'-aft:len0') then
        i = i + 1
        call get_command_argument(i,tag,status=ios)
        if (ios.ne.0) then
            write(0,*) 'minage: error reading argument for -aft:len0 LEN1,LEN2,...'
            write(0,*) 'Received error flag: ios=',ios
            write(0,*) 'Did you remember to include the list of initial fission track lengths?'
            call error_exit(1)
        endif
        aft_n0 = 1
        do
            j = index(tag,',')
            if (j.ne.0) then
                tag(j:j) = ' '
                aft_n0 = aft_n0 + 1
            else
                exit
            endif
        enddo
        allocate(aft_len0(aft_n0))
        read(tag,*,iostat=ios) (aft_len0(j),j=1,aft_n0)
        if (ios.ne.0) then
            write(0,*) 'minage: error reading fission track lengths with -aft:len0 LEN1,LEN2,...'
            write(0,*) 'Tried to parse "',trim(tag),'" as lengths'
            call error_exit(1)
        endif

    elseif (trim(tag).eq.'-aft:segmentation') then
        i = i + 1
        call get_command_argument(i,tag,status=ios)
        if (tag.eq.'ON'.or.tag.eq.'on'.or.tag.eq.'1') then
            doSegmentationCarlson1990 = .true.
        elseif (tag.eq.'OFF'.or.tag.eq.'off'.or.tag.eq.'0') then
            doSegmentationCarlson1990 = .false.
        else
            call usage('minage: error reading segmentation ON/OFF flag')
        endif

    elseif (trim(tag).eq.'-aft:userbias') then
        i = i + 1
        call get_command_argument(i,tag,status=ios)
        if (tag.eq.'ON'.or.tag.eq.'on'.or.tag.eq.'1') then
            doEtchingUserBiasCorrectionWillett1997 = .true.
        elseif (tag.eq.'OFF'.or.tag.eq.'off'.or.tag.eq.'0') then
            doEtchingUserBiasCorrectionWillett1997 = .false.
        else
            call usage('minage: error reading etching/user bias ON/OFF flag')
        endif


    !*********************************************!
    ! Apatite (U-Th)/He age calculation variables !
    !*********************************************!
    elseif (trim(tag).eq.'-ahe:radius') then
        i = i + 1
        call get_command_argument(i,tag,status=ios)
        if (ios.ne.0) then
            write(0,*) 'minage: error reading argument for -ahe:radius R1,R2,...'
            write(0,*) 'Received error flag: ios=',ios
            write(0,*) 'Did you remember to include the list of grain radii?'
            call error_exit(1)
        endif
        ahe_nradius = 1
        do
            j = index(tag,',')
            if (j.ne.0) then
                tag(j:j) = ' '
                ahe_nradius = ahe_nradius + 1
            else
                exit
            endif
        enddo
        allocate(ahe_radius(ahe_nradius))
        read(tag,*,iostat=ios) (ahe_radius(j),j=1,ahe_nradius)
        if (ios.ne.0) then
            write(0,*) 'minage: error reading apatite radii with -ahe:radius R1,R2,...'
            write(0,*) 'Tried to parse "',trim(tag),'" as radii'
            call error_exit(1)
        endif

    elseif (trim(tag).eq.'-ahe:verbosity') then
        i = i + 1
        call get_command_argument(i,tag,status=ios)
        read(tag,*) ahe_verbosity
    elseif (trim(tag).eq.'-ahe:fd') then
        call usage('*****FD*****')

    elseif (trim(tag).eq.'-ahe:dtfactor') then              ! Timestep resampling
        i = i + 1
        call get_command_argument(i,tag,status=ios)
        read(tag,*) ahe_dt_var
    elseif (trim(tag).eq.'-ahe:dtma') then
        i = i + 1
        call get_command_argument(i,tag,status=ios)
        read(tag,*) ahe_dt_var
        ahe_dt_var = -abs(ahe_dt_var)
    elseif (trim(tag).eq.'-ahe:nnodes') then                ! Number of spatial nodes
        i = i + 1
        call get_command_argument(i,tag,status=ios)
        read(tag,*) ahe_nnodes
    elseif (trim(tag).eq.'-ahe:beta') then                  ! Implicitness weight
        i = i + 1
        call get_command_argument(i,tag,status=ios)
        read(tag,*) ahe_beta
    elseif (trim(tag).eq.'-ahe:taumax') then                ! Max dimensionless time to retain He
        i = i + 1
        call get_command_argument(i,tag,status=ios)
        read(tag,*) ahe_taumax



    elseif (trim(tag).eq.'-advanced') then
        call usage('*****ADVANCED*****')

    else
        call usage('minage: no option '//trim(tag))
    endif

    if (ios.ne.0) then
        call usage('minage: error parsing "'//trim(tag)//'" flag arguments')
    endif

    i = i + 1

enddo

return
end subroutine



!--------------------------------------------------------------------------------------------------!

subroutine usage(str)
implicit none
character(len=*) :: str
integer :: i, j
logical :: printAdvancedOptions
logical :: printFiniteDifference
printAdvancedOptions = .false.
printFiniteDifference = .false.
i = index(str,'*****ADVANCED*****')
if (i.ne.0) then
    printAdvancedOptions = .true.
endif
j = index(str,'*****FD*****')
! write(0,*) 'j',j
if (j.ne.0) then
    printAdvancedOptions = .true.
    printFiniteDifference = .true.
    i = j
endif
if (str.ne.'') then
    if (printAdvancedOptions) then
        write(0,*) trim(str(1:i-1))
    else
        write(0,*) trim(str)
    endif
    write(0,*)
endif
write(0,*) 'Usage: minage -temp READTQTEC_TEMP_FILE [-dep READTQTEC_DEP_FILE] [...options...]'
write(0,*)
if (printAdvancedOptions) then
write(0,*) '******************** BASIC MINAGE OPTIONS *******************'
endif
write(0,*) '-temp READTQTEC_TEMP_FILE  Temperature history output from readtqtec'
write(0,*) '-dep READTQTEC_DEP_FILE    Depth history (default: track age from beginning of model)'
write(0,*) '-aft AFT_FILE              Apatite fission track age'
write(0,*) '-aft:history DT FILE       Apatite fission track distribution every DT Ma'
write(0,*) '-ahe AHE_FILE              Apatite (U-Th)/He age'
write(0,*) '-advanced                  See advanced options'
write(0,*)
if (printAdvancedOptions) then
write(0,*)
write(0,*) '******************** ADVANCED AFT OPTIONS *******************'
write(0,*) '-aft:len0 LEN1,LEN2,...    Lengths of fission tracks generated every timestep (microns)'
write(0,*) '    [5 15-um, 7 16-um, 7 17-um, 1 18-um => 20 FTs, avg len = 16.2 um]'
write(0,*) '-aft:segmentation ON|OFF   Turn Carlson (1990) track segmentation [on]/off'
write(0,*) '-aft:userbias ON|OFF       Turn Willett (1997) user bias/etching correction [on]/off'
write(0,*)
write(0,*)
write(0,*) '******************** ADVANCED AHE OPTIONS *******************'
write(0,*) '-ahe:radius R1,R2,...  Apatite grain radii to calculate ages [20,40,...,120,140]'
write(0,*) '-ahe:verbosity LVL     Verbosity level for apatite helium calculation [0]'
write(0,*) '-ahe:dtfactor VAL      Timestep resampling factor [0.1]'
write(0,*) '    MinAge will try to reduce dt_resamp to 0.5*dr^2/max_diffusivity, but will use the'
write(0,*) '    maximum of that and the following values:'
write(0,*) '        VAL > 0: dt_resamp = dt_init * VAL'
write(0,*) '        VAL < 0: dt_resamp = |VAL|'
write(0,*) '-ahe:dtma DTMA         Resampled timestep (alternative to -ahe:dtfactor, sets VAL = -|DTMA|)'
write(0,*) '-ahe:nnodes NNODES     NNODES = Number of spatial nodes + 2 BC nodes [102]'
write(0,*) '-ahe:beta BETA         Finite difference implicitness coefficient [0.85]'
write(0,*) '-ahe:taumax TAUMAX     Maximum dimensionless time to retain He [0.4]'
write(0,*) '-ahe:fd                Print description of finite difference options'
write(0,*)
endif
if (printFiniteDifference) then
write(0,*) '******************** FINITE DIFFERENCE PARAMETERS ********************'
write(0,*) '* The default finite difference parameters (listed below) are selected to produce:'
write(0,*) '*     1. An accurate solution to the diffusion-production equation over a wide range of'
write(0,*) '*        thermochronological scenarios'
write(0,*) '*     2. A solution that remains stable in those situations'
write(0,*) '*     3. A calculation completes "relatively quickly" (within a few seconds)'
write(0,*) '*'
write(0,*) '* Depending on your temperature history and finite difference parameters, you may still'
write(0,*) '* experience long runtimes (not too common) or stability issues (more common).'
write(0,*) '* In general, increasing the number of nodes (reducing node spacing) improves the'
write(0,*) '* solution accuracy at the expense of stability and runtime. Decreasing the timestep'
write(0,*) '* size increases solution accuracy and stability at the expense of runtime. Stability'
write(0,*) '* also depends on the node spacing relative to the rate of helium diffusion, which'
write(0,*) '* is a function of temperature.'
write(0,*) '*'
write(0,*) '* MinAge provides several options to adjust the finite difference setup that will'
write(0,*) '* affect the accuracy, speed, and stability of the calculation:'
write(0,*) '*     1. Resampled timestep size (-ahe:dtfactor or -ahe:dtma)'
write(0,*) '*     2. Number of nodes (-ahe:nnodes)'
write(0,*) '*     3. Implicitness weight (-ahe:beta)'
write(0,*) '*'
write(0,*) '* Solution stability typically becomes a problem at high temperatures, when the amount'
write(0,*) '* of helium moving from node to node is very high. The most accurate way to handle this'
write(0,*) '* problem is to decrease the timestep size, but at high temperatures a stable timestep'
write(0,*) '* can lead to prohibitively long runtimes. An alternative simplifying assumption is to'
write(0,*) '* set a threshold diffusion timescale above which all helium is lost (set to zero).'
write(0,*) '* In MinAge, we calculate the dimensionless time, tau, for 1 Ma:'
write(0,*) '*     tau = diffusivity*(1 Ma)/radius^2'
write(0,*) '* And compare this to the value of tau_max set on the command line:'
write(0,*) '*     4. Threshold dimensionless time for complete helium escape (-ahe:taumax)'
write(0,*) '*'
write(0,*) '* Default Parameter Values'
write(0,*) '*   Resampled timestep size:  VAL=0.1'
write(0,*) '*                             dt_resamp=max(0.5*dr^2/max_diffusivity,VAL*dt_input)'
write(0,*) '*   Number of spatial nodes:  102'
write(0,*) '*   Implicitness coefficient: 0.85'
write(0,*) '*   Retention threshold time: 0.4 (dimensionless)'
write(0,*)
endif
call error_exit(1)
stop
end subroutine
