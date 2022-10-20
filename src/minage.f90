module minage

character(len=512) :: readtqtec_temp_file

end module


!==================================================================================================!


program main

use minage, only: readtqtec_temp_file

use fission_track, only: calc_fission_track_distribution, &
                         calc_fission_track_ages, &
                         ft_nbins, &
                         ft_hist_corr, &
                         ft_age_corr, &
                         ft_retention_age_corr

implicit none

double precision :: time_ma
double precision :: dt_ma
double precision, allocatable :: temp_array(:,:)
integer :: i
integer :: ios
integer :: n
integer :: ihorizon, nhorizons
integer, parameter :: max_nhorizons = 100
real :: dum(max_nhorizons)
character(len=512) :: line
double precision, allocatable :: ages_ft_length(:)
double precision, allocatable :: ages_ft_count(:)
double precision, allocatable :: temp_ft_retention_age_corr(:)
integer, allocatable :: hists_ft(:,:)
character(len=32) :: fmt_string



!call read_temp_history() ! PUT THE STUFF AFTER THIS INTO A SUBROUTINE

! Read and parse temperature history from readtqtec -temp output
readtqtec_temp_file = 'MWX_test6.3.temp'
open(unit=21,file=readtqtec_temp_file,status='old')

! Total time, sample rate, and number of temperature-time points
read(21,*,iostat=ios) time_ma, dt_ma, n
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
do i = 1,max_nhorizons
    read(line,*,iostat=ios) dum(1:i)
    if (ios.ne.0) then
        nhorizons = i-1
        exit
    endif
enddo
backspace(21)

! Allocate and initialize the temperature array
allocate(temp_array(nhorizons,n))
temp_array = 0.0d0

! Fill the temperature array
do i = 1,n
    read(21,*,iostat=ios) (temp_array(ihorizon,i),ihorizon=1,nhorizons)
    if (ios.ne.0) then
        write(0,*) 'minage: error reading line',n,' of tqtec output temperature file'
        call error_exit(1)
    endif
enddo



!----
! Fission track distributions and ages
!----
! Allocate memory to fission track arrays
allocate(ages_ft_length(nhorizons))
allocate(ages_ft_count(nhorizons))
allocate(temp_ft_retention_age_corr(nhorizons))
ages_ft_length = 0.0d0
ages_ft_count = 0.0d0

! Calculate fission track distributions and ages for each horizon
do ihorizon = 1,nhorizons

    ! Calculate fission track distributions (variables in fission_track module)
    !     ft_nbins               number of bins in the fission track distributions
    !     ft_dlen                width of histogram bins
    !     ft_hist_raw            uncorrected for track segmentation
    !     ft_hist                corrected for track segmentation
    !     ft_hist_corr           corrected for etching and user bias
    ! print *,'calculating fission track distribution'
    call calc_fission_track_distribution(temp_array(ihorizon,:), n, dt_ma)
    if (.not.allocated(hists_ft)) then
        allocate(hists_ft(nhorizons,ft_nbins))
    endif
    hists_ft(ihorizon,1:ft_nbins) = ft_hist_corr(1:ft_nbins)

    ! Calculate fission track ages (variables in fission_track module)
    !     ft_age_corr            age based on fission track lengths, corrected (see module for details)
    !     ft_retention_age_corr  age based on number of tracks, corrected (see module for details)
    ! print *,'calculating fission track ages'
    call calc_fission_track_ages(n,dt_ma)
    ages_ft_length(ihorizon) = ft_age_corr
    ages_ft_count(ihorizon) = ft_retention_age_corr

    ! Save temperature at time of retention age
    i = nint(-time_ma/dt_ma - ft_retention_age_corr/dt_ma)
    temp_ft_retention_age_corr(ihorizon) = temp_array(ihorizon,i)
enddo
close(21)

write(*,*) 'CORRECTED FT AGE'
write(fmt_string,'("("I6,"F8.3)")') nhorizons
write(*,fmt_string) (ages_ft_length(i),i=1,nhorizons)
write(*,*) 'CORRECTED RETENTION AGE'
write(*,fmt_string) (ages_ft_count(i),i=1,nhorizons)
write(*,*) 'TEMP AT TIME OF CORRECTED RETENTION AGE'
write(fmt_string,'("("I6,"F8.2)")') nhorizons
write(*,fmt_string) (temp_ft_retention_age_corr(i),i=1,nhorizons)
write(*,*) 'FISSION TRACK LENGTH HISTOGRAM (CORRECTED FOR SEGMENTATION AND ETCHING/USER BIAS)'
write(fmt_string,'("(I6,"I6,"I8)")') nhorizons
do i = 1,ft_nbins
     write(*,fmt_string) i,(hists_ft(ihorizon,i),ihorizon=1,nhorizons)
enddo

end program