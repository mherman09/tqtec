module fission_track
!----
! Apatite fission track modeling module
!
! Originally written by Matt Legg at Penn State
! Adapted into Modern Fortran by Matt Herman at CSU Bakersfield
!----

! Physical constants
double precision, parameter :: boltzmann_constant = 3.297623483d-27 ! [kcal/K]             K
double precision, parameter :: univ_gas_constant = 1.98720425864d-3 ! [kcal/mol/K]         R
double precision, parameter :: planck_constant = 1.583d-37          ! [kcal*s]             H

! Exponent in power-law for initial radial defect distribution
double precision, parameter :: defect_power_law_constant = 0.206d0  !                      N

! Activation energy to eliminate defect
double precision, parameter :: defect_activation_energy = 40.6d0    ! [kcal/mol]           Q

! Rate constant for axial shortening
double precision, parameter :: axial_shortening_constant = 1.81d0   ! [micron]             A

! Average initial fission track length (SHOULD THIS BE A VARIABLE?)
double precision, parameter :: ft_len_avg_init = 16.2d0             ! [micron?]            LOAVE

! Ratio of present-day average length for spontaneous tracks in Durango apatite to LOAVE (Donelick and Miller, 1991)
double precision, parameter :: PST = 0.893d0                        !                      PST


! Fission track distribution parameters (from calc_fission_track_distribution)
integer :: ft_nbins
double precision :: ft_hist_len_max
double precision :: ft_hist_dlen
double precision, allocatable :: ft_hist(:)
double precision, allocatable :: ft_hist_corr(:)



!==================================================================================================!
contains
!==================================================================================================!


subroutine calc_fission_track_length(len_init, time_ma, temp_celsius, n, dt_ma, len_final)
!----
! Determine fission track lengths over a temperature-time history, based on Carlson (1990)
!
! Inputs:
!   len_init: initial track length (microns)
!   time_ma: n-point double precision array of times (Ma)
!   temp_celsius: n-point double precision array of temperatures (C)
!   dt_ma: time step size (Ma)
!
! Outputs:
!   len_final: n-point double precision array of track lengths at each time (microns)
!----

implicit none

! Arguments
double precision :: len_init
integer :: n
double precision :: time_ma(n)
double precision :: temp_celsius(n)
double precision :: dt_ma
double precision :: len_final(n)

! Local variables
double precision :: time_seconds(n)
double precision :: temp_kelvin(n)
double precision :: dt_seconds
double precision :: annealing(n)
double precision :: annealing_sum(n)
double precision :: arg
double precision :: sum
double precision :: const
integer :: i


! Convert time and temperature units
time_seconds = time_ma*1d6*365.25d0*24.0d0*60.0d0*60.0d0
temp_kelvin = temp_celsius + 273.0d0
dt_seconds = dt_ma*1d6*365.25d0*24.0d0*60.0d0*60.0d0


!----
! Carlson (1990) Equation 4 describes length of a fission track from a given thermal history:
!
!                            n                                        n
!                      ( k )      [ St (           (  -Q    )      )]
!     l_as = l_0 - A * (---)    * [ S  ( T(y) * exp(--------) * dy )]
!                      ( h )      [ S0 (           ( R*T(y) )      )]
!
!         l_as: etchable length after axial shortening
!         l_0: initial, unannealed length
!         A: Empirical rate constant for axial shortening
!         k: Boltzmann's constant
!         h: Planck's constant
!         n: Exponent in power law expression for initial radial defect distribution
!         t: Time of interest
!         T: Temperature
!         Q: Activation energy for atomic motions that eliminate defects
!         R: Universal gas constant
!         y: Dummy time variable for integration

! Annealing for each time step (integrand function)
do i = 1,n
    arg = -defect_activation_energy/(univ_gas_constant*temp_kelvin(i))
    annealing(i) = temp_kelvin(i) * dt_seconds * exp(arg)
enddo

! Calculate cumulative sum of annealing terms from n to 1 (estimating integral from 0 to t)
sum = annealing(n)
annealing_sum(n) = sum
do i = n-1,1,-1
    sum = annealing_sum(i+1) + annealing(i)
    annealing_sum(i) = sum
enddo

! Calculate the final length (finishing the formula)
const = axial_shortening_constant * (boltzmann_constant/planck_constant)**defect_power_law_constant
do i = 1,n
    len_final(i) = len_init - const * annealing_sum(i)**defect_power_law_constant

    ! Set negative lengths to zero
    if (len_final(i).lt.0.0d0) then
        len_final(i) = 0.0d0
    endif
enddo

return
end subroutine


!--------------------------------------------------------------------------------------------------!


subroutine calc_fission_track_distribution(time_ma, temp_celsius, n, dt_ma)
!----
! Determine fission track length distribution for a temperature-time history
!
! Inputs:
!   time_ma: n-point double precision array of times (Ma)
!   temp_celsius: n-point double precision array of temperatures (C)
!   dt_ma: time step size (Ma)
!
! Outputs (fission_track module variables)
!   ft_hist_len_max: maximum length included in fission track length histogram
!   ft_hist_dlen:    fission track length histogram bin width
!   ft_nbins:        number of bins in fission track length histogram
!   ft_hist:         fission track length histogram
!   ft_hist_corr:    corrected fission track length histogram

! References:
!     Green et al. (1986) - Initial track length distribution
!     Carlson (1990) - Kinetics of track annealing
!     Willett (1997) - Correction for etching/user bias
!----

implicit none

! Arguments
integer :: n
double precision :: time_ma(n)
double precision :: temp_celsius(n)
double precision :: dt_ma

! Local variables
integer, parameter :: ninit = 20
double precision :: len_init(ninit)
double precision :: len_final(ninit,n)
double precision :: len
double precision :: ratio
integer :: ibin
integer :: i, j


! Set initial fission track distribution (Green et al., 1986)
len_init(1:5) = 15.0d0
len_init(6:12) = 16.0d0
len_init(13:19) = 17.0d0
len_init(20) = 18.0d0


! Calculate fission track lengths for the input temperature history
do i = 1,ninit
    call calc_fission_track_length(len_init(i), time_ma, temp_celsius, n, dt_ma, len_final(i,:))
enddo


!----
! Generate track length histogram
!----
! Initialize histogram dimensions
ft_hist_len_max = 20.0d0                                            ! maximum track length for histogram (microns)
ft_hist_dlen = 1.0d0                                                ! histogram bin width (microns)
ft_nbins = int((ft_hist_len_max+0.5d0*ft_hist_dlen)/ft_hist_dlen)   ! number of bins

! Allocate memory to histogram arrays and initialize bin counts to 0
allocate(ft_hist(ft_nbins))
allocate(ft_hist_corr(ft_nbins))
ft_hist = 0
ft_hist_corr = 0

! Load raw track lengths into histogram array
do i = 1,ninit
    do j = 1,n
        ! Bin array index
        ibin = int((len_final(i,j)+0.5d0*ft_hist_dlen)/ft_hist_dlen)
        ! Update count
        if (0.lt.ibin .and. ibin.le.ft_nbins) then
            ft_hist(ibin) = ft_hist(ibin) + 1
        endif
    enddo
enddo

! Empirical count adjustments to distributions based on track lengths
do i = 1,ft_nbins

    ! Track length for the bin
    len = dble(i)*ft_hist_dlen

    ! Account for segmentation of lengths up to 11 microns
    ! Carlson (1990) empirical fit shown in Figure 8
    if (len.le.11.0d0) then
        ft_hist(i) = ft_hist(i) - int((-0.1d0*len + 1.2d0)*ft_hist(i))
    endif
    if (ft_hist(i).lt.0) then
        ft_hist(i) = 0
    endif

    ! Correct histogram for etching/user bias based on Willett (1997) Equation 4
    if (len.le.6.0d0) then
        ft_hist_corr(i) = 0
    elseif (len.le.11.0d0) then
        ratio = len/16.2d0          ! ratio of bin length to average initial length (16.2 microns)
        ft_hist_corr(i) = ft_hist(i) * ((2.862d0*ratio)-1.2104d0)
    endif
    if (ft_hist_corr(i).lt.0) then
        ft_hist_corr(i) = 0
    endif
enddo

return
end subroutine


!--------------------------------------------------------------------------------------------------!


subroutine calc_fission_track_ages(hist,hist_corr,nbins)
implicit none

integer :: nbins
integer :: hist(nbins)
integer :: hist_corr(nbins)

integer :: ntracks, ntracks_corr
integer :: i


! Check whether histograms have been calculated
if (.not.allocated(ft_hist)) then
    write(0,*) 'calc_fission_track_ages: fission track distributions must be calculated first'
endif

! Calculate sum of track lengths from computed uncorrected and corrected distributions
ntracks = 0
ntracks_corr = 0
do i = 1,nbins
    ntracks = ntracks + hist(i)
    ntracks_corr = ntracks_corr + hist(i)
enddo

! Calculate sum of track lengths multiplied by time step

return
end subroutine

end module