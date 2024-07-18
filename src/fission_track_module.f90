!--------------------------------------------------------------------------------------------------!
! Module Fission Track                                                                             !
!                                                                                                  !
! Authors:                                                                                         !
!     - Matt Herman (current Modern Fortran version, i.e., what you are looking at right now!)     !
!     - Matt Legg (original version; PSU MS thesis)                                                !
!                                                                                                  !
! This module contains variables and subroutines for the generation of fission tracks produced by  !
! spontaneous fission in geological materials and the analysis of ages corresponding to fission    !
! track analysis.                                                                                  !
!                                                                                                  !
! References                                                                                       !
! Carlson, W. (1990). Mechanisms and kinetics of apatite fission-track annealing. American         !
!     Mineralogist, 75, 1120–1139.                                                                 !
! Donelick, R.A., Miller, D.S. (1991). Enhanced tint fission track densities in low spontaneous    !
!     track density apatites using 252Cf-derived fission fragment tracks: A model and experimental !
!     observations. International Journal of Radiation Applications and Instrumentation. Part D.   !
!     Nuclear Tracks and Radiation Measurements, 18, 3, 301-307.                                   !
! Green, P.F., Duddy, I.R., Gleadow, A.J.W., Tingate, P.R., Laslett, G.M. (1986). Thermal          !
!     annealing of fission tracks in apatite: 1. A qualitative description. Chemical Geology:      !
!     Isotope Geoscience Section, 59, 237-253.                                                     !
! Ketcham, R. A. (2005). Forward and Inverse Modeling of Low-Temperature Thermochronometry Data.   !
!     Reviews in Mineralogy and Geochemistry, 58(1), 275–314.                                      !
! Legg, M.J. (2010). The Tectonic and Thermal Evolution of Hawke's Bay Basin, New Zealand. Penn    !
!     State MS Thesis. 72 pp.                                                                      !
! Willett, S.D. (1997). Inverse modeling of annealing of fission tracks in apatite: 1. A           !
!     controlled random search method. American Journal of Science, 297, 10, 939-969.              !
!--------------------------------------------------------------------------------------------------!



module fission_track


    ! Exponent in power-law for initial radial defect distribution
    double precision, parameter :: defect_power_law_constant = 0.206d0


    ! Activation energy to eliminate defect
    double precision, parameter :: defect_activation_energy = 40.6d0    ! [kcal/mol]


    ! Rate constant for axial shortening
    double precision, parameter :: axial_shortening_constant = 1.81d0   ! [micron]


    ! Average initial fission track length
    double precision, parameter :: ft_len_avg_init = 16.2d0             ! [micron]


    ! Ratio of present-day average length for spontaneous tracks in Durango apatite to the average
    ! initial fission track length (Donelick and Miller, 1991)
    double precision, parameter :: PST = 0.893d0



    ! Fission track distribution parameters
    integer, parameter :: ft_ninit = 20                                 ! number of FTs per timestep
    double precision, allocatable :: ft_len_final_all(:,:)              ! FT lengths at each timestep
    integer :: ft_nbins                                                 ! # FT bins in distribution
    double precision :: ft_hist_len_max                                 ! max FT length in distribution
    double precision :: ft_hist_dlen                                    ! bin width in distribution
    integer, allocatable :: ft_hist(:)
    integer, allocatable :: ft_hist_raw(:)
    integer, allocatable :: ft_hist_corr(:)



    ! Fission track ages (from calc_fission_track_ages)
    double precision :: ft_age                         ! AGE3
    double precision :: ft_age_corr                    ! FTAGE
    double precision :: ft_retention_age_raw           ! AGE0
    double precision :: ft_retention_age               ! AGE1
    double precision :: ft_retention_age_corr          ! AGE2




contains !-----------------------------------------------------------------------------------------!




    subroutine generate_ft_len_temp_history(nft_init, &
                                            len_init, &
                                            nt, &
                                            temp_celsius, &
                                            dep_km, &
                                            dt_ma, &
                                            ft_len)

    !----
    ! Calculate fission track lengths due to axial shortening of tracks produced over a
    ! temperature-time history. Fission tracks of length <len_init> are generated every timestep,
    ! and Equation 4 from Carlson (1990) is used to determine the final length of each track.
    !
    ! Carlson (1990) Equation 4 describes length of a fission track from a given thermal history
    ! as it undergoes axial shortening:
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
    !
    ! Inputs:
    !   nft_init:     number of fission tracks generated each fission event
    !   len_init:     nft_init-point double precision initial fission track length array (microns)
    !   nt:           number of time steps
    !   temp_celsius: n-point double precision temperature history array (Celsius)
    !   dep_km:       n-point double precision depth history array (km)
    !   dt_ma:        double precision timestep size (Ma)
    !
    ! Outputs:
    !   ft_len:       (nft_init x nt) double precision fission track length array (microns)
    !----


    use physical_constants, only: boltzmann_constant, &
                                planck_constant, &
                                universal_gas_constant


    implicit none


    ! Arguments
    integer :: nft_init
    double precision :: len_init(nft_init)
    integer :: nt
    double precision :: temp_celsius(nt)
    double precision :: dep_km(nt)
    double precision :: dt_ma
    double precision :: ft_len(nft_init,nt)


    ! Local variables
    integer :: i
    integer :: j
    double precision :: dt_seconds
    double precision :: temp_kelvin(nt)
    double precision :: annealing(nt)
    double precision :: annealing_sum(nt)
    double precision :: arg
    double precision :: sum
    double precision :: const
    logical :: keepTracks



    ! Convert time and temperature units
    temp_kelvin = temp_celsius + 273.0d0
    dt_seconds = dt_ma * 1d6 * 365d0 * 24.0d0 * 60.0d0 * 60.0d0 ! SHOULD YEAR BE 365.25?



    ! Annealing for each time step (integrand function)
    do i = 1,nt
        arg = -defect_activation_energy / (universal_gas_constant * temp_kelvin(i))
        annealing(i) = temp_kelvin(i) * dt_seconds * exp(arg)
    enddo



    ! Calculate cumulative sum of annealing terms from n to 1 (estimating integral from 0 to t)
    ! Going in reverse order simulates a fission track produced every timestep that goes through
    ! temperature history from that time until end of temperature history.
    sum = annealing(nt)
    annealing_sum(nt) = sum
    do i = nt-1,1,-1
        sum = annealing_sum(i+1) + annealing(i)
        annealing_sum(i) = sum
    enddo



    ! Calculate the present-day lengths of fission tracks generated at each timestep
    const = axial_shortening_constant * (boltzmann_constant/planck_constant)**defect_power_law_constant
    do j = 1,nft_init
        do i = 1,nt

        ! Finish Carlson (1990) Equation 4
            ft_len(j,i) = len_init(j) - const * annealing_sum(i)**defect_power_law_constant

        ! Set negative lengths to zero
            if (ft_len(j,i).lt.0.0d0) then
                ft_len(j,i) = 0.0d0
            endif

        enddo
    enddo



    ! Remove tracks from time period before mineral actually formed
    keepTracks = .false. ! At start of run, assume mineral has not formed and do not keep tracks
    do i = 1,nt          ! For each time...
        if (.not.keepTracks) then
            if (dep_km(i).le.0.0d0) then
                keepTracks = .true.     ! Horizon below surface, mineral exists, keep tracks from here on
            else
                ft_len(:,i) = 0.0d0     ! Horizon above surface, remove tracks corresponding to this time
            endif
        endif
    enddo



    return

    end subroutine generate_ft_len_temp_history



!--------------------------------------------------------------------------------------------------!



    subroutine segment_ft_distribution(nbins, binwid, len_min, hist, hist_corr)

    !----
    ! Add effects of fission track segmentation to a track length distribution. Segmentation occurs
    ! due to radial shortening of tracks ~12 microns and shorter (Carlson, 1990). Each track splits
    ! into two shorter tracks of random length. Figure 8 from Carlson (1990) shows the number of
    ! segmented tracks as a function of axial length, with fit given by:
    !
    !     f = s * (len_sg - len_as)
    !
    !     len_sg:    segmentation initiation length (empirically 12 microns)
    !     len_as:    mean annealed length from axial shortening alone
    !     s:         empirically derived slope from fraction segmented vs. length plot (Fig. 8)
    !
    ! Inputs:
    !   nbins:      number of histogram bins
    !   binwid:     histogram bin width
    !   len_min:    minimum length in histogram
    !   hist:       input histogram array
    !
    ! Outputs:
    !   hist_corr:  histogram with segmented tracks
    !----


    implicit none


    ! Arguments
    integer, intent(in) :: nbins
    double precision, intent(in) :: binwid
    double precision, intent(in) :: len_min
    integer, intent(in) :: hist(nbins)
    integer, intent(out) :: hist_corr(nbins)


    ! Local variables
    integer :: i
    integer :: j
    integer :: ct
    double precision :: len
    double precision :: fraction_segmented
    double precision :: arg
    double precision :: ran
    double precision, parameter :: s = 0.10d0
    integer, allocatable :: seed(:)
    integer :: time(3)
    character(len=8) :: string
    double precision :: x
    integer :: ibin



    ! Initialize a random number generator (for determining lengths of segmented tracks)
    if (.not.allocated(seed)) then
        allocate(seed(1))
        call itime(time)                                          ! itime() returns hr, mn, sc to int array
        write(string,'(I0.2I0.2I0.2)') time(1), time(2), time(3)  ! combine to single integer
        read(string,*) seed(1)
        call random_seed(PUT=seed)
    endif



    ! For each length bin...
    do i = 1,nbins

        ! Fission track length in the bin
        len = dble(i-1)*binwid + len_min


        ! Calculate fraction of tracks that are segmented at this length (Carlson, 1990)
        if (len.le.12.0d0) then
            fraction_segmented = s * (12.0d0 - len)
        else
            fraction_segmented = 0.0d0
            hist_corr(i) = hist(i)
            cycle
        endif


        ! For each track in this bin, determine whether it is segmented or not
        ct = 0
        do j = 1,hist(i)

            arg = dble(j)*fraction_segmented - dble(ct)   ! test for segmenting this track


            if (arg.ge.1.0d0) then

                !***** SEGMENT TRACK *****!

                call random_number(ran) !-----------------!---- WHERE TO SEGMENT?

                x = len*ran !-----------------------------!---- SEGMENT 1
                ibin = int((x+0.5d0*binwid)/binwid)       ! which bin is segmented track 1 in?
                if (ibin.lt.1) then                       !
                    ibin = 1                              !
                elseif (ibin.gt.nbins) then               !
                    ibin = nbins                          !
                endif                                     !
                hist_corr(ibin) = hist_corr(ibin) + 1     ! add segmented track 1 to bin

                x = len - x !-----------------------------!---- SEGMENT 2
                ibin = int((x+0.5d0*binwid)/binwid)       ! which bin is segmented track 1 in?
                if (ibin.lt.1) then                       !
                    ibin = 1                              !
                elseif (ibin.gt.nbins) then               !
                    ibin = nbins                          !
                endif                                     !
                hist_corr(ibin) = hist_corr(ibin) + 1     ! add segmented track 1 to bin

                hist_corr(i) = hist_corr(i) - 1 !---------!---- REMOVE INITIAL TRACK
            endif

            ct = ct + 1

        enddo


        ! Set negative counts to zero
        if (hist_corr(i).lt.0) then
            hist_corr(i) = 0
        endif


    enddo

    return

    end subroutine segment_ft_distribution



!--------------------------------------------------------------------------------------------------!



    subroutine correct_ft_distribution_etching_userbias(nbins, &
                                                        binwid, &
                                                        len_min, &
                                                        hist, &
                                                        len0, &
                                                        hist_corr)

    !----
    ! Correct histogram for preferential etching and user bias counting of longer fission tracks.
    ! Green (1988) and Willett (1997) showed that fission track density (on which the fission track
    ! age is based) is not 1:1 linearly proportional to mean fission track length. Rather, density
    ! of longer tracks is proportional to their lengths, but shorter tracks go to zero density.
    ! Equation 4 of Willett (1997) relates the track density to the track lengths:
    !
    !             len/len0 <= 0.423:   p/p0 = 0
    !     0.423 < len/len0 <= 0.650:   p/p0 = 2.862*len/len0 - 1.2104
    !     0.650 < len/len0:            p/p0 = len/len0
    !
    !         len:  track length (microns)
    !         len0: initial mean track length (microns)
    !         p:    correected track density
    !         p0:   control sample track density
    !
    ! Inputs:
    !   nbins:      number of histogram bins
    !   binwid:     histogram bin width (microns)
    !   len_min:    minimum fission track length in the histogram (microns)
    !   hist:       fission track length histogram
    !   len0:       initial average fission track length (microns)
    !
    ! Outputs
    !   hist_corr:  corrected fission track length histogram
    !----


    implicit none


    ! Arguments
    integer, intent(in) :: nbins
    double precision, intent(in) :: binwid
    double precision, intent(in) :: len_min
    integer, intent(in) :: hist(nbins)
    double precision, intent(in) :: len0
    integer, intent(out) :: hist_corr(nbins)


    ! Local variables
    integer :: i
    double precision :: len
    double precision :: ratio
    double precision :: a
    double precision :: b



    ! Domain of piecewise function
    a = 0.423d0*len0
    b = 0.650d0*len0


   ! Correct count in each bin of the histogram
    do i = 1,nbins

        ! Track length for the bin
        len = dble(i-1)*binwid + len_min


        ! Correct histogram for etching and user bias based on Willett (1997) Equation 4
        if (len.le.a) then
            ! No short tracks
            hist_corr(i) = 0
        elseif (len.le.b) then
            ! Reduce number of medium length tracks
            ratio = len/len0   ! ratio of bin length to average initial length
            hist_corr(i) = int( dble(hist(i)) * ((2.862d0*ratio)-1.2104d0) )
        else
            ! Keep number of long tracks
            hist_corr(i) = hist(i)
        endif


        ! Set negative counts to zero
        if (hist_corr(i).lt.0) then
            hist_corr(i) = 0
        endif

    enddo


    return

    end subroutine correct_ft_distribution_etching_userbias





!--------------------------------------------------------------------------------------------------!



    subroutine calc_fission_track_ages(n,dt_ma)
    !----
    ! Determine fission track ages
    !
    ! Inputs (arguments)
    !   n:                       number of time steps
    !   dt_ma:                   time step size (Ma)
    !
    ! Inputs (fission_track module variables)
    !   ft_len_final_all:        array of fission track lengths
    !   ft_hist:                 fission track length histogram
    !   ft_hist_corr:            corrected fission track length histogram
    !   ft_nbins:                number of bins in fission track length histogram
    !   ft_hist_dlen:            fission track length histogram bin width
    !
    ! Outputs (fission_track module variables)
    !   ft_retention_age_raw:    track-count age before any corrections
    !   ft_retention_age:        track-count age after segmentation correction
    !   ft_retention_age_corr:   track-count age after segmentation + etching/user bias corrections
    !   ft_age:                  track-length age before any corrections
    !   ft_age_corr:             time step size (Ma)
    !----


    implicit none


    ! Arguments
    integer :: n
    double precision :: dt_ma


    ! Local variables
    integer :: ntracks_raw
    integer :: ntracks
    integer :: ntracks_corr
    double precision :: total_track_len
    double precision :: total_track_len_hist
    double precision :: total_track_len_hist_corr
    double precision :: len
    double precision :: ratio
    integer :: i, j



    ! Check whether histograms have been calculated
    if (.not.allocated(ft_hist).or..not.allocated(ft_hist_corr)) then
        write(0,*) 'calc_fission_track_ages: fission track distributions must be calculated first'
        write(0,*) 'with calc_fission_track_distribution()'
        call error_exit(1)
    endif


    !********** RETENTION AGES **********!
    ! Retention ages are based on the number of tracks counted:
    !   - Total number of tracks produced over history: total_time / dt_ma * ft_ninit
    !   - Number of tracks remaining at present: ntracks
    !   - Apparent age: ntracks * dt_ma / ft_init    <-- RETENTION AGE

    ! Calculate total number of tracks in uncorrected and corrected distributions
    ntracks_raw = 0                                                         ! TOT0 (not considering segmentation)
    ntracks = 0                                                             ! TOT1 (segmentation)
    ntracks_corr = 0                                                        ! TOT2 (segmentation + etching/user bias)
    do i = 1,ft_nbins
        ntracks_raw  = ntracks_raw  + ft_hist_raw(i)
        ntracks      = ntracks      + ft_hist(i)
        ntracks_corr = ntracks_corr + ft_hist_corr(i)
    enddo


    ! Calculate retention ages (= #tracks * timestep / total #tracks generated)
    ft_retention_age_raw  = ntracks_raw  * dt_ma/ft_ninit                   ! AGE0
    ft_retention_age      = ntracks      * dt_ma/ft_ninit                   ! AGE1
    ft_retention_age_corr = ntracks_corr * dt_ma/ft_ninit                   ! AGE2



    !********** LENGTH DISTRIBUTION AGES *********!
    ! Length distribution ("fission track") ages are based on the track length distributions:
    !   - Total length of tracks produced over history: total_time / dt_ma * ft_len_avg_init * ft_ninit
    !   - Total length of tracks remaining at present: total_track_len
    !   - Apparent age: total_track_len * dt_ma / ft_len_avg_init / ft_ninit    <-- FISSION TRACK AGE

    ! Calculate total sum of all track lengths produced in model
    total_track_len = 0.0d0                                                 ! TOT3 (intermediate step)
    do i = 1,ft_ninit
        do j = 1,n
            total_track_len = total_track_len + ft_len_final_all(i,j)
        enddo
    enddo


    ! Divide by average track length and number of tracks; multiply by time step duration
    total_track_len = total_track_len*dt_ma/dble(ft_ninit)/ft_len_avg_init  ! TOT3 (final step)


    ! Calculate raw (no corrections) fission track age (Ketcham, 2005; Equation 14)
    ft_age = total_track_len/PST                                            ! AGE3


    ! Total track lengths for corrected fission track distributions
    total_track_len_hist = 0.0d0                                            ! FTOT5 (segmentation only)
    total_track_len_hist_corr = 0.0d0                                       ! FTOT2 (segmentation + etching/user bias)


    ! Calculate total fission track lengths for corrected distributions
    do i = 1,ft_nbins
        len = dble(i)*ft_hist_dlen
        total_track_len_hist      = total_track_len_hist      + ft_hist(i)*len
        total_track_len_hist_corr = total_track_len_hist_corr + ft_hist_corr(i)*len
    enddo


    ! Ratio between segmentation + etching/user bias corrected and segmentation corrected total lengths
    ratio = total_track_len_hist_corr/total_track_len_hist


    ! Correct fission track age by multiplying raw age by ratio
    ft_age_corr = ft_age * ratio                                            ! FTAGE


    return
    end subroutine




end module