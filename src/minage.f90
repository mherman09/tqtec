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
character(len=512) :: aft_file
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
integer :: aft_n0
double precision, allocatable :: aft_len0(:)
logical :: doSegmentationCarlson1990
logical :: doEtchingUserBiasCorrectionWillett1997

! Apatite (U-Th)/He variables
double precision :: ahe_beta
double precision :: ahe_dt_var
integer :: ahe_nnodes
double precision :: ahe_taumax



end module


!==================================================================================================!


program main


use minage, only: readtqtec_temp_file, &
                  readtqtec_dep_file, &
                  aft_file, &
                  ahe_file, &
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
if (aft_file.ne.'') then
    call calc_aft_ages()
endif


! Apatite (U-Th)/He (variables and descriptions in apatite_helium_module.f90)
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



!----
!   Generate fission tracks and calculate ages
!----
do ihorizon = 1,nhorizons ! for each horizon...


    ! Calculate fission track lengths for the temperature history -> ftlen
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


    ! Load the fission track lengths into a histogram
    call load_histogram(aft_n0*ntimes, aft_len, binwid, xmin, nbins, aft_hist(:,ihorizon))
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



! Print results to AFT file
open(unit=22,file=aft_file,status='unknown')

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


! Free array memory
deallocate(aft_len)
deallocate(aft_hist)
deallocate(aft_hist_corr)
deallocate(aft_age)
deallocate(aft_retention_age)


return

end subroutine




!--------------------------------------------------------------------------------------------------!




subroutine calc_ahe_ages()
!----
! Calculate apatite (U-Th)/He ages from a temperature-depth history. Subroutines and variables
! for apatite-helium calculations can be found in apatite_helium_module.f90
!----


use minage, only: ahe_file, &
                  nhorizons, &
                  dt_ma, &
                  ntimes, &
                  temp_celsius_array, &
                  dep_km_array, &
                  ahe_beta, &
                  ahe_dt_var, &
                  ahe_nnodes, &
                  ahe_taumax

use apatite_helium, only: calc_apatite_he_age


implicit none


integer :: ihorizon                                              !
integer :: iradius                                               !
integer, parameter :: nradius = 7                                !
double precision :: grain_radius_array(nradius)                  !
double precision :: radius_microns                               !
double precision :: he_age                                       !
double precision, allocatable :: horizon_ahe_age(:,:)            !
character(len=32) :: fmt_string                                  !



write(*,*) 'minage: calculating apatite (U-Th)/He ages'


! Open apatite-helium output file
open(unit=23,file=ahe_file,status='unknown')



! Apatite grain radii
grain_radius_array(1) = 20.0d0
grain_radius_array(2) = 40.0d0
grain_radius_array(3) = 60.0d0
grain_radius_array(4) = 80.0d0
grain_radius_array(5) = 100.0d0
grain_radius_array(6) = 120.0d0
grain_radius_array(7) = 140.0d0



! Array for grain-size-based ages in each tracked horizon
if (.not.allocated(horizon_ahe_age)) then
    allocate(horizon_ahe_age(nhorizons,nradius))
endif
horizon_ahe_age = 0.0d0


! For each horizon...
do ihorizon = 1,nhorizons


    ! For each grain radius...
    do iradius = 1,nradius


        ! if (printProgress) then
        !     write(*,2301) 'minage: working on horizon',ihorizon,'of',nhorizons, &
        !                 ': grain size',iradius,'of',nradius
        !     2301 format(X,2(A,I6,X,A,I6))
        ! endif


        ! Calculate (U-Th)/He age for this temperature-depth-time history and grain radius
        ! See apatite_helium_module.f90 for correction details
        radius_microns = grain_radius_array(iradius)
        call calc_apatite_he_age(temp_celsius_array(ihorizon,:), &
                                 dep_km_array(ihorizon,:), &
                                 ntimes, &
                                 dt_ma, &
                                 ahe_dt_var, &
                                 radius_microns, &
                                 ahe_nnodes, &
                                 ahe_beta, &
                                 ahe_taumax, &
                                 he_age)


        ! Save grain-size-based ages for this horizon
        horizon_ahe_age(ihorizon,iradius) = he_age

    enddo

enddo


! Print results to apatite helium file
write(23,*) '# Apatite (U-Th)/He ages (Ma)'
write(23,*) '# Radius(microns) Age(Ma)'
write(fmt_string,'("(F12.1,"I6,"F12.3)")') nhorizons
do iradius = 1,nradius
    write(23,fmt_string) grain_radius_array(iradius), &
                         (horizon_ahe_age(ihorizon,iradius),ihorizon=1,nhorizons)
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

use minage, only: readtqtec_temp_file, &
                  readtqtec_dep_file, &
                  aft_file, &
                  ahe_file, &
                  isOutputDefined, &
                  aft_n0, &
                  aft_len0, &
                  doSegmentationCarlson1990, &
                  doEtchingUserBiasCorrectionWillett1997, &
                  ahe_beta, &
                  ahe_dt_var, &
                  ahe_nnodes, &
                  ahe_taumax

implicit none

integer :: i
integer :: j
integer :: narg
integer :: ios
character(len=512) :: tag


! Initialize variables
readtqtec_temp_file = ''
readtqtec_dep_file = ''
aft_file = ''
ahe_file = ''
isOutputDefined = .false.


! Apatite fission track variables
aft_n0 = 0
doSegmentationCarlson1990 = .true.
doEtchingUserBiasCorrectionWillett1997 = .true.

! Apatite (U-Th)/He variables
ahe_beta = 0.85d0
ahe_dt_var = 0.1d0
ahe_nnodes = 102
ahe_taumax = 0.40d0


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


    ! Output files
    elseif (trim(tag).eq.'-aft') then
        i = i + 1
        call get_command_argument(i,aft_file,status=ios)
        isOutputDefined = .true.

    elseif (trim(tag).eq.'-ahe') then
        i = i + 1
        call get_command_argument(i,ahe_file,status=ios)
        isOutputDefined = .true.


    ! Apatite fission track age calculation variables
    elseif (trim(tag).eq.'-aft:len0') then
        i = i + 1
        call get_command_argument(i,tag,status=ios)
        if (ios.ne.0) then
            call usage('minage: error reading argument for -aft:len0 LEN1,LEN2,...')
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
            call usage('minage: error reading initial fission track length for -aft:len0')
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


    ! Apatite (U-Th)/He age calculation variables
    elseif (trim(tag).eq.'-ahe:beta') then
        i = i + 1
        call get_command_argument(i,tag,status=ios)
        read(tag,*) ahe_beta
    elseif (trim(tag).eq.'-ahe:dtfactor') then
        i = i + 1
        call get_command_argument(i,tag,status=ios)
        read(tag,*) ahe_dt_var
    elseif (trim(tag).eq.'-ahe:dtma') then
        i = i + 1
        call get_command_argument(i,tag,status=ios)
        read(tag,*) ahe_dt_var
        ahe_dt_var = -ahe_dt_var
    elseif (trim(tag).eq.'-ahe:nnodes') then
        i = i + 1
        call get_command_argument(i,tag,status=ios)
        read(tag,*) ahe_nnodes
    elseif (trim(tag).eq.'-ahe:taumax') then
        i = i + 1
        call get_command_argument(i,tag,status=ios)
        read(tag,*) ahe_taumax



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
integer :: i
logical :: printAdvancedOptions
i = index(str,'*****ADVANCED*****')
if (i.ne.0) then
    printAdvancedOptions = .true.
    str(i:i+18) = ' '
else
    printAdvancedOptions = .false.
endif
if (str.ne.'') then
    write(0,*) trim(str)
    write(0,*)
endif
write(0,*) 'Usage: minage -temp READTQTEC_TEMP_FILE [-dep READTQTEC_DEP_FILE] [...options...]'
write(0,*)
write(0,*) '-temp READTQTEC_TEMP_FILE  Temperature history output from readtqtec'
write(0,*) '-dep READTQTEC_DEP_FILE    Depth history (default: track age from beginning of model)'
write(0,*) '-aft AFT_FILE              Apatite fission track age'
write(0,*) '-ahe AHE_FILE              Apatite (U-Th)/He age'
write(0,*) '-advanced'
write(0,*)
write(0,*) 'Apatite Fission Track Options:'
write(0,*) '-aft:len0 LEN1,LEN2,...    Lengths of fission tracks generated every timestep (microns)'
write(0,*) '-aft:segmentation ON|OFF   Turn Carlson (1990) track segmentation correction [on]/off'
write(0,*) '-aft:userbias ON|OFF       Turn Willett (1997) user bias/etching correction [on]/off'
write(0,*)
write(0,*) 'Apatite (U-Th)/He Options:'
write(0,*) '-ahe:radius R1,R2,...  Apatite grain radii to calculate ages (20,40,...,120,140)'
write(0,*) '-ahe:param             Print finite difference parameter values'
write(0,*)
write(0,*)
if (printAdvancedOptions) then
write(0,*) '******************** ADVANCED OPTIONS *******************'
write(0,*)
write(0,*) 'FINITE DIFFERENCE'
write(0,*) 'The default finite difference parameters were selected to produce (1) an accurate'
write(0,*) 'solution to the diffusion-production equation over a wide range of geological scenarios'
write(0,*) 'that (2) remains stable, and (3) completes relatively quickly. Increasing the number of'
write(0,*) 'nodes improves accuracy, at the expense of stability and runtime. Stability depends on'
write(0,*) 'the node spacing relative to the rate of helium diffusion; we provide several options'
write(0,*) 'below to adjust the calculation. This can be done by modifying the details of the finite'
write(0,*) 'difference setup:'
write(0,*) '    1. Decrease the timestep size'
write(0,*) '    2. Decrease the number of nodes'
write(0,*) '    3. Increase the weighting of the implicitness in the calculation'
write(0,*) 'Alternatively, the assumptions for diffusion in the geological system can be modified:'
write(0,*) '    4. Define a threshold temperature for complete escape of the diffusion product'
write(0,*)
write(0,*) 'Default Parameter Values'
write(0,*) '  Number of spatial nodes:  102 (including 2 boundary condition nodes)'
write(0,*) '  Resampled timestep size:  0.5*(node_spacing)^2/max_diffusivity (minval=0.1*dt_input)'
write(0,*) '  Implicitness coefficient: 0.85'
write(0,*) '  Retention threshold:      diffusivity*(1 Ma)/grain_radius'
write(0,*)
write(0,*) '-ahe:nnodes NNODES     Number of spatial nodes + 2 BC nodes'
write(0,*) '-ahe:dtma DTMA         Resampled timestep size, in Ma'
write(0,*) '-ahe:beta BETA         Finite difference implicitness coefficient'
write(0,*) '-ahe:temp TEMP         Maximum temperature to retain He'
write(0,*)
endif
call error_exit(1)
stop
end subroutine
