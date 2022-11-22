module minage

character(len=512) :: readtqtec_temp_file
character(len=512) :: aft_file
character(len=512) :: ahe_file
logical :: isOutputDefined

integer :: nhorizons
integer, parameter :: nhorizons_max = 100

double precision :: time_ma
double precision :: dt_ma
integer :: ntimes
double precision, allocatable :: temp_celsius_array(:,:)

! Apatite (U-Th)/He variables
double precision :: ahe_beta
double precision :: ahe_dt_max_reduction_factor
integer :: ahe_nnodes
double precision :: ahe_taumax

end module


!==================================================================================================!


program main

use minage, only: readtqtec_temp_file, &
                  aft_file, &
                  ahe_file, &
                  isOutputDefined !, &
                !   nhorizons, &
                !   time_ma, &
                !   dt_ma, &
                !   ntimes, &
                !   temp_celsius_array

! use apatite_helium, only: calc_apatite_he_age
! use radiogenic_helium, only: atomic_mass_th232, &
!                              atomic_mass_u235, &
!                              atomic_mass_u238, &
!                              decay_th232, &
!                              decay_u235, &
!                              decay_u238, &
!                              calc_he_production_rate, calc_u_th_he_age


implicit none

! double precision :: time_ma
! double precision :: dt_ma
! double precision, allocatable :: temp_array(:,:)
! integer :: i
! integer :: ios
! integer :: n
! double precision, allocatable :: ages_ft_length(:)
! double precision, allocatable :: ages_ft_count(:)
! double precision, allocatable :: temp_ft_retention_age_corr(:)
! integer, allocatable :: hists_ft(:,:)
! double precision :: radius_microns

! double precision :: volume_sphere, mass_u238, mass_u235, mass_th232, production_rate, mol_he4, duration, age
! double precision :: mol_u238, mol_u235, mol_th232, t, dt


!----
! Parse command line arguments
!----
call gcmdln()


!----
! Input error checks
!----
! Is input temperature history defined?
if (readtqtec_temp_file.eq.'') then
    call usage('minage: no readtqtec output temperature file defined')
endif

! Is an output defined?
if (.not.isOutputDefined) then
    call usage('minage: no output defined')
endif


!----
! Read temperature history
!----
call read_temp_history()



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

subroutine calc_aft_ages()

use minage, only: aft_file, &
                  nhorizons, &
                  time_ma, &
                  dt_ma, &
                  ntimes, &
                  temp_celsius_array

use fission_track, only: calc_fission_track_distribution, &
                         calc_fission_track_ages, &
                         ft_hist_dlen, &
                         ft_nbins, &
                         ft_hist_corr, &
                         ft_age_corr, &
                         ft_retention_age_corr


implicit none

integer :: ihorizon
integer :: itemp
integer :: ibin
integer, allocatable :: horizon_ft_hist_corr(:,:)
double precision, allocatable :: horizon_ft_length_age(:)
double precision, allocatable :: horizon_ft_count_age(:)
double precision, allocatable :: temp_ft_retention_age_corr(:)
character(len=32) :: fmt_string


write(*,*) 'minage: calculating apatite fission track ages'
open(unit=22,file=aft_file,status='unknown')


! For each horizon...
do ihorizon = 1,nhorizons

    ! Calculate fission track distribution
    call calc_fission_track_distribution(temp_celsius_array(ihorizon,:), ntimes, dt_ma)

    ! Save corrected fission track distribution for this horizon
    if (.not.allocated(horizon_ft_hist_corr)) then
        allocate(horizon_ft_hist_corr(nhorizons,ft_nbins))
    endif
    horizon_ft_hist_corr(ihorizon,1:ft_nbins) = ft_hist_corr(1:ft_nbins)

    ! Calculate fission track ages
    call calc_fission_track_ages(ntimes,dt_ma)

    ! Save length-based and count-based ages for this horizon
    if (.not.allocated(horizon_ft_length_age)) then
        allocate(horizon_ft_length_age(nhorizons))
    endif
    if (.not.allocated(horizon_ft_count_age)) then
        allocate(horizon_ft_count_age(nhorizons))
    endif
    horizon_ft_length_age(ihorizon) = ft_age_corr
    horizon_ft_count_age(ihorizon) = ft_retention_age_corr

    ! Save temperature at time of count-based (retention) age
    if (.not.allocated(temp_ft_retention_age_corr)) then
        allocate(temp_ft_retention_age_corr(nhorizons))
    endif
    itemp = nint(-time_ma/dt_ma - ft_retention_age_corr/dt_ma)
    temp_ft_retention_age_corr(ihorizon) = temp_celsius_array(ihorizon,itemp)
enddo


! Print results to AFT file
write(22,*) 'Corrected track-length-based apatite fission track ages (Ma)'
write(fmt_string,'("(6X,"I6,"F8.3)")') nhorizons
write(22,fmt_string) (horizon_ft_length_age(ihorizon),ihorizon=1,nhorizons)
write(22,*)

write(22,*) 'Corrected track-count-based apatite fission track ages (Ma)'
write(22,fmt_string) (horizon_ft_count_age(ihorizon),ihorizon=1,nhorizons)
write(22,*)

write(22,*) 'Temperatures (C) at time of track-count-based (retention) ages'
write(fmt_string,'("(6X,"I6,"F8.2)")') nhorizons
write(22,fmt_string) (temp_ft_retention_age_corr(ihorizon),ihorizon=1,nhorizons)
write(22,*)

write(22,*) 'Fission track length histograms (corrected for segmentation & etching/user bias)'
write(fmt_string,'("(F6.1,"I6,"I8)")') nhorizons
do ibin = 1,ft_nbins
    write(22,fmt_string) ibin*ft_hist_dlen, &
                         (horizon_ft_hist_corr(ihorizon,ibin),ihorizon=1,nhorizons)
enddo


! Close file
close(22)


return

end subroutine


!--------------------------------------------------------------------------------------------------!


subroutine calc_ahe_ages()

use minage, only: ahe_file, &
                  nhorizons, &
                  dt_ma, &
                  ntimes, &
                  temp_celsius_array, &
                  ahe_beta, &
                  ahe_dt_max_reduction_factor, &
                  ahe_nnodes

use apatite_helium, only: calc_apatite_he_age

implicit none

integer :: ihorizon
integer :: iradius
integer, parameter :: nradius = 1
double precision :: grain_radius_array(nradius)
double precision :: radius_microns
double precision :: he_age
double precision, allocatable :: horizon_ahe_age(:,:)
character(len=32) :: fmt_string


! Apatite grain radii to calculate diffusion
! grain_radius_array(1) = 10.0d0
! grain_radius_array(2) = 20.0d0
! grain_radius_array(3) = 50.0d0
! grain_radius_array(4) = 100.0d0
! grain_radius_array(5) = 200.0d0
! grain_radius_array(6) = 500.0d0
! grain_radius_array(7) = 1000.0d0
grain_radius_array(1) = 100.0d0

write(*,*) 'minage: calculating apatite (U-Th)/He ages'
open(unit=23,file=ahe_file,status='unknown')


! Allocate array to save grain-size-based ages for this horizon
if (.not.allocated(horizon_ahe_age)) then
    allocate(horizon_ahe_age(nhorizons,nradius))
endif
horizon_ahe_age = 0.0d0


! For each horizon...
do ihorizon = 2,2

    ! For each grain radius...
    do iradius = 1,1

        write(*,*) 'minage: working on horizon',ihorizon,' and grain size',iradius,' of',nradius

        ! Calculate (U-Th)/He age
        radius_microns = grain_radius_array(iradius)
        call calc_apatite_he_age(temp_celsius_array(ihorizon,:), &
                                 ntimes, &
                                 dt_ma, &
                                 ahe_dt_max_reduction_factor, &
                                 radius_microns, &
                                 ahe_nnodes, &
                                 ahe_beta, &
                                 he_age)

        ! Save grain-size-based ages for this horizon
        horizon_ahe_age(ihorizon,iradius) = he_age
    enddo

enddo


! Print results to AFT file
write(23,*) 'Apatite (U-Th)/He ages (Ma)'
write(fmt_string,'("(F6.1,"I6,"F8.3)")') nhorizons
do iradius = 1,nradius
    write(23,fmt_string) grain_radius_array(iradius), &
                         (horizon_ahe_age(ihorizon,iradius),ihorizon=1,nhorizons)
enddo
! write(22,fmt_string) (horizon_ft_length_age(ihorizon),ihorizon=1,nhorizons)
! write(22,*)

! write(22,*) 'Corrected track-count-based apatite fission track ages (Ma)'
! write(22,fmt_string) (horizon_ft_count_age(ihorizon),ihorizon=1,nhorizons)
! write(22,*)

! write(22,*) 'Temperatures (C) at time of track-count-based (retention) ages'
! write(fmt_string,'("(6X,"I6,"F8.2)")') nhorizons
! write(22,fmt_string) (temp_ft_retention_age_corr(ihorizon),ihorizon=1,nhorizons)
! write(22,*)

! write(22,*) 'Fission track length histograms (corrected for segmentation & etching/user bias)'


! Close file
close(23)


return

end subroutine

!--------------------------------------------------------------------------------------------------!


subroutine gcmdln()
!----
! Parse command line
!----

use minage, only: readtqtec_temp_file, &
                  aft_file, &
                  ahe_file, &
                  isOutputDefined, &
                  ahe_beta, &
                  ahe_dt_max_reduction_factor, &
                  ahe_nnodes, &
                  ahe_taumax

implicit none

integer :: i
integer :: narg
integer :: ios
character(len=512) :: tag


! Initialize variables
readtqtec_temp_file = ''
aft_file = ''
ahe_file = ''
isOutputDefined = .false.

! Apatite (U-Th)/He variables
ahe_beta = 0.85d0
ahe_dt_max_reduction_factor = 0.1d0
ahe_nnodes = 502


! Count number of command line arguments
narg = command_argument_count()
if (narg.eq.0) then
    call usage('')
endif


! Read and parse command line
i = 1
do while (i.le.narg)

    call get_command_argument(i,tag)

    if (trim(tag).eq.'-temp') then
        i = i + 1
        call get_command_argument(i,readtqtec_temp_file,status=ios)

    elseif (trim(tag).eq.'-aft') then
        i = i + 1
        call get_command_argument(i,aft_file,status=ios)
        isOutputDefined = .true.

    elseif (trim(tag).eq.'-ahe') then
        i = i + 1
        call get_command_argument(i,ahe_file,status=ios)
        isOutputDefined = .true.


    ! Apatite (U-Th)/He variables
    elseif (trim(tag).eq.'-ahe:beta') then
        i = i + 1
        call get_command_argument(i,tag,status=ios)
        read(tag,*) ahe_beta
    elseif (trim(tag).eq.'-ahe:dtfactor') then
        i = i + 1
        call get_command_argument(i,tag,status=ios)
        read(tag,*) ahe_dt_max_reduction_factor
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
write(0,*) 'Usage: minage -temp TQTEC_TEMP_FILE [...options...]'
write(0,*)
write(0,*) '-temp TQTEC_TEMP_FILE  Temperature history output from readtqtec'
write(0,*) '-aft AFT_FILE          Apatite fission track age'
write(0,*) '-ahe AHE_FILE          Apatite (U-Th)/He age'
write(0,*)
write(0,*) 'Apatite (U-Th)/He options:'
write(0,*) '-ahe:beta BETA         Finite difference implicitness coefficient (0.85)'
write(0,*) '-ahe:nnodes NNODES     Number of spatial nodes + 2 BC nodes (502)'
write(0,*) '-ahe:taumax TAUMAX     Maximum dimensionless time to retain He for 1 Ma (0.30)'
write(0,*) '-ahe:dtfactor FACTOR   Timestep resampling factor (0.10)'
write(0,*)
call error_exit(1)
stop
end subroutine
