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


use minage, only: aft_file, &
                  nhorizons, &
                  time_ma, &
                  dt_ma, &
                  ntimes, &
                  temp_celsius_array, &
                  dep_km_array


use fission_track, only:       &
    generate_fts_carlson_1990, &
    segment_fts_carlson_1990,  &
    correct_fts_etching_userbias_willett_1997, &
                         load_histogram, &
                         calc_ft_retention_age, &
                         calc_ft_age
!                          ft_hist_dlen, &
!                          ft_nbins, &
!                          ft_hist_corr, &
!                          ft_age_corr, &
!                          ft_retention_age_corr



implicit none


! Local variables
integer :: i
integer :: ihorizon                                              ! Horizon index
integer :: nft0
double precision, allocatable :: ftlen0(:)
double precision :: len0
double precision, allocatable :: ftlen(:,:)
double precision :: xmin
double precision :: binwid
integer :: ibin
integer :: nbins
integer, allocatable :: fthist(:)
double precision :: ft_age
double precision :: ft_retention_age
! integer :: itemp                                                 ! Temperature index
character(len=32) :: fmt_string
integer, allocatable :: fthist_all(:,:)
! double precision, allocatable :: horizon_ft_length_age(:)        ! FT-length-based age
! double precision, allocatable :: horizon_ft_count_age(:)         ! FT-count-based age
! double precision, allocatable :: temp_ft_retention_age_corr(:)   ! 



write(*,*) 'minage: calculating apatite fission track ages'


! Open fission track output file
open(unit=22,file=aft_file,status='unknown')




! Lengths of fission tracks generated during a fission event
nft0 = 20
allocate(ftlen0(nft0))
allocate(ftlen(nft0,ntimes))
ftlen0(1:5) = 15.0d0
ftlen0(6:12) = 16.0d0
ftlen0(13:19) = 17.0d0
ftlen0(20) = 18.0d0
len0 = 0.0d0
do i = 1,nft0
    len0 = len0 + ftlen0(i)
enddo
len0 = len0/dble(nft0)
write(0,*) 'len0 = ',len0

! Histogram parameters
xmin = 1.0d0
binwid = 1.0d0
nbins = 20
allocate(fthist(nbins))


! For each horizon...
do ihorizon = 1,nhorizons


    ! Calculate fission track lengths for the temperature history
    call generate_fts_carlson_1990(nft0, &
                                      ftlen0, &
                                      ntimes, &
                                      temp_celsius_array(ihorizon,:), &
                                      dep_km_array(ihorizon,:), &
                                      dt_ma, &
                                      'green-et-al-1986', &
                                      ftlen)

    ! Load the fission track lengths into a histogram
    call load_histogram(nft0*ntimes, ftlen, binwid, xmin, nbins, fthist)

    ! Calculate fission track ages
    call calc_ft_retention_age(nbins, fthist, nft0, dt_ma, ft_retention_age)
    write(0,*) 'ft_retention_age = ',ft_retention_age
    call calc_ft_age(nbins, binwid, xmin, fthist, nft0, len0, dt_ma, ft_age)
    write(0,*) 'ft_age = ',ft_age
    
!     ! Save corrected fission track distribution for this tracked horizon
!     ! See fission_track_module.f90 for correction details
    if (.not.allocated(fthist_all)) then
        allocate(fthist_all(nhorizons,nbins))
    endif
    fthist_all(ihorizon,1:nbins) = fthist(1:nbins)


!     ! Calculate a series of fission track ages:
!     ! See subroutines in fission_track_module.f90 for different age calculations
!     call calc_fission_track_ages(ntimes,dt_ma)


!     ! Save length-based and count-based ages for this horizon
!     if (.not.allocated(horizon_ft_length_age)) then
!         allocate(horizon_ft_length_age(nhorizons))
!     endif
!     if (.not.allocated(horizon_ft_count_age)) then
!         allocate(horizon_ft_count_age(nhorizons))
!     endif
!     horizon_ft_length_age(ihorizon) = ft_age_corr                ! Length-distribution-based age
!     horizon_ft_count_age(ihorizon) = ft_retention_age_corr       ! Track-count-based age


!     ! Save temperature at time of track-count-based (retention) age
!     if (.not.allocated(temp_ft_retention_age_corr)) then
!         allocate(temp_ft_retention_age_corr(nhorizons))
!     endif
!     if (time_ma.lt.0.0d0) then
!         itemp = nint(-time_ma/dt_ma - ft_retention_age_corr/dt_ma)
!     else
!         itemp = nint(time_ma/dt_ma - ft_retention_age_corr/dt_ma)
!     endif
!     if (itemp.lt.1) then
!         itemp = 1
!     elseif (itemp.gt.ntimes) then
!         itemp = ntimes
!     endif
!     temp_ft_retention_age_corr(ihorizon) = temp_celsius_array(ihorizon,itemp)


enddo   ! end of loop over horizons


! ! Print results to AFT file
! write(22,*) 'Corrected track-length-based apatite fission track ages (Ma)'
! write(fmt_string,'("(6X,"I6,"F8.3)")') nhorizons
! write(22,fmt_string) (horizon_ft_length_age(ihorizon),ihorizon=1,nhorizons)
! write(22,*)


! write(22,*) 'Corrected track-count-based apatite fission track ages (Ma)'
! write(22,fmt_string) (horizon_ft_count_age(ihorizon),ihorizon=1,nhorizons)
! write(22,*)


! write(22,*) 'Temperatures (C) at time of track-count-based (retention) ages'
! write(fmt_string,'("(6X,"I6,"F8.2)")') nhorizons
! write(22,fmt_string) (temp_ft_retention_age_corr(ihorizon),ihorizon=1,nhorizons)
! write(22,*)


write(22,*) 'Fission track length histograms'
write(fmt_string,'("(F6.1,"I6,"I8)")') nhorizons
do ibin = 1,nbins
    write(22,fmt_string) xmin+dble(ibin-1)*binwid, &
                         (fthist_all(ihorizon,ibin),ihorizon=1,nhorizons)
enddo




! Close file
close(22)


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
                  ahe_beta, &
                  ahe_dt_var, &
                  ahe_nnodes, &
                  ahe_taumax

implicit none

integer :: i
integer :: narg
integer :: ios
character(len=512) :: tag


! Initialize variables
readtqtec_temp_file = ''
readtqtec_dep_file = ''
aft_file = ''
ahe_file = ''
isOutputDefined = .false.


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
write(0,*) '-advanced                  See all advanced options'
write(0,*)
write(0,*) 'Apatite (U-Th)/He Options:'
write(0,*) '-ahe:radius R1,R2,...  Apatite grain radii to calculate ages (20,40,...,120,140)'
write(0,*) '-ahe:param             Print finite difference parameter values'
write(0,*)
write(0,*)
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
call error_exit(1)
stop
end subroutine
