!--------------------------------------------------------------------------------------------------!
! Module Fission Track                                                                             !
!                                                                                                  !
! Authors:                                                                                         !
!     - Matt Herman (current Modern Fortran version, i.e., what you are looking at right now!)     !
!     - Matt Legg (original version; PSU MS thesis)                                                !
!     - Kevin Furlong (original version; supervisor)                                               !
!                                                                                                  !
! This module contains variables and subroutines for the generation of fission tracks produced by  !
! spontaneous fission in geological materials and the analysis of ages corresponding to fission    !
! track analysis.                                                                                  !
!                                                                                                  !
! References                                                                                       !
! Carlson, W. (1990). Mechanisms and kinetics of apatite fission-track annealing. American         !
!     Mineralogist, 75, 1120–1139.                                                                 !
! Crowley, K.D. (1993). LENMODEL: A forward model for calculating length distributions and         !
!     fission-track ages in apatite. Computers & Geosciences 19, 5, 619-626.                       !
! Donelick, R.A., Miller, D.S. (1991). Enhanced tint fission track densities in low spontaneous    !
!     track density apatites using 252Cf-derived fission fragment tracks: A model and experimental !
!     observations. International Journal of Radiation Applications and Instrumentation. Part D.   !
!     Nuclear Tracks and Radiation Measurements, 18, 3, 301-307.                                   !
! Green, P.F. (1988). The relationship between track shortening and fission track age reduction in !
!     apatite: Combined influences of inherent instability, annealing anisotropy, length bias and  !
!     system calibration. Earth and Planetary Science Letters, 89(3–4), 335–352.                   !
! Green, P.F., Duddy, I.R., Gleadow, A.J.W., Tingate, P.R., Laslett, G.M. (1986). Thermal          !
!     annealing of fission tracks in apatite: 1. A qualitative description. Chemical Geology:      !
!     Isotope Geoscience Section, 59, 237-253.                                                     !
! Ketcham, R. A. (2005). Forward and Inverse Modeling of Low-Temperature Thermochronometry Data.   !
!     Reviews in Mineralogy and Geochemistry, 58(1), 275–314.                                      !
! Ketcham, R.A., Donelick, R.A., Carlson, W.D. (1999). Variability of apatite fission-track        !
!     annealing kinetics: III. Extrapolation to geological time scales. American Mineralogist, 84, !
!     1235-1255.                                                                                   !
! Legg, M.J. (2010). The Tectonic and Thermal Evolution of Hawke's Bay Basin, New Zealand.         !
!     Pennsylvania State University MS Thesis. 72 pp.                                              !
! Willett, S.D. (1997). Inverse modeling of annealing of fission tracks in apatite: 1. A           !
!     controlled random search method. American Journal of Science, 297, 10, 939-969.              !
!--------------------------------------------------------------------------------------------------!



module fission_track


    ! Fission track generation subroutines
    PUBLIC :: generate_fts_carlson_1990
    PUBLIC :: generate_fts_ketcham_et_al_1999
    PUBLIC :: segment_fts_carlson_1990
    PUBLIC :: correct_fts_etching_userbias_willett_1997

    ! Fission track age subroutines
    PUBLIC :: calc_ft_age
    PUBLIC :: calc_ft_retention_age
    PRIVATE :: total_len




contains



!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!-------------------------------- FISSION TRACK ANNEALING ROUTINES --------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!



    subroutine generate_fts_carlson_1990( &
        nft_init,                         &
        len_init,                         &
        nt,                               &
        temp_celsius,                     &
        dep_km,                           &
        dt_ma,                            &
        kinetic_par,                      &
        ft_len                            &
    )

    !----
    ! Calculate fission track lengths due to axial shortening of tracks produced over a
    ! temperature-time history. <nft_init> fission tracks of length <len_init(nft_init)> are
    ! generated every timestep, and Equation 4 from Carlson (1990) is used to determine the final
    ! length of each track.
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
    !   kinetic_par:  string indicating kinetic parameters
    !
    ! Outputs:
    !   ft_len:       (nft_init x nt) double precision fission track length array (microns)
    !----


    use physical_constants, only: boltzmann_constant, &
                                  planck_constant, &
                                  universal_gas_constant


    implicit none


    ! Arguments
    integer, intent(in) :: nft_init
    double precision, intent(in) :: len_init(nft_init)
    integer, intent(in) :: nt
    double precision, intent(in) :: temp_celsius(nt)
    double precision, intent(in) :: dep_km(nt)
    double precision, intent(in) :: dt_ma
    character(len=*), intent(in) :: kinetic_par
    double precision, intent(out) :: ft_len(nft_init,nt)


    ! Local variables
    integer :: i
    integer :: j
    double precision :: defect_power_law_constant
    double precision :: defect_activation_energy
    double precision :: axial_shortening_constant
    double precision :: dt_seconds
    double precision :: temp_kelvin(nt)
    double precision :: annealing(nt)
    double precision :: annealing_sum(nt)
    double precision :: arg
    double precision :: sum
    double precision :: const
    logical :: keepTracks



    ! Define kinetic parameters for annealing (Carlson, 1990, Table 2)
    if (kinetic_par.eq.'green-et-al-1986') then
        defect_power_law_constant = 0.206d0
        defect_activation_energy = 40.6d0   ! [kcal/mol]
        axial_shortening_constant = 1.81d0  ! [micron]
    elseif (kinetic_par.eq.'donelick-1988') then
        defect_power_law_constant = 0.111d0
        defect_activation_energy = 49.0d0   ! [kcal/mol]
        axial_shortening_constant = 5.85d0  ! [micron]
    elseif (kinetic_par.eq.'carlson-1990-composite') then
        defect_power_law_constant = 0.141d0
        defect_activation_energy = 46.5d0   ! [kcal/mol]
        axial_shortening_constant = 4.66d0  ! [micron]
    else
        write(0,*) 'generate_ft_len_temp_history: error setting kinetic parameters for annealing'
        write(0,*) 'Read kinetic_par =',trim(kinetic_par)
        write(0,*) 'Options are:'
        write(0,*) '    "green-et-al-1986"'
        write(0,*) '    "donelick-1988"'
        write(0,*) '    "carlson-1990-composite"'
        call error_exit(1)
    endif


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
    do i = 1,nt
        if (.not.keepTracks) then
            if (dep_km(i).le.0.0d0) then
                keepTracks = .true.     ! Below surface, mineral exists, keep tracks from here on
            else
                ft_len(:,i) = 0.0d0     ! Above surface, remove tracks from this time
            endif
        endif
    enddo



    return

    end subroutine generate_fts_carlson_1990



!--------------------------------------------------------------------------------------------------!



    subroutine generate_fts_ketcham_et_al_1999( &
        nft_init,                               &
        len_init,                               &
        nt,                                     &
        temp_celsius,                           &
        dep_km,                                 &
        dt_ma,                                  &
        kinetic_par,                            &
        ft_len                                  &
    )

    !----
    ! Calculate fission track lengths due to axial shortening of tracks produced over a
    ! temperature-time history. <nft_init> fission tracks of length <len_init(nft_init)> are
    ! generated every timestep, and the "equivalent time" backwards solution method of Crowley
    ! (1993) is used with the parameters derived in Ketcham et al. (1999) to determine the final
    ! length of each track.
    !
    ! The equation to solve (Ketcham et al., 1999) is either the "fanning Arrhenius model":
    !
    !     ( 1 - r^beta )^alpha
    !     ( ---------- )
    !     (    beta    )                     (  ln(t) - C2   )
    !     --------------------   = C0 + C1 * (-------------- )
    !            alpha                       (  (1/T) - C3   )
    !
    ! or the "curvilinear Arrhenius model":
    !
    !     ( 1 - r^beta )^alpha
    !     ( ---------- )
    !     (    beta    )                     (  ln(t) - C2  )
    !     --------------------   = C0 + C1 * (------------- )
    !            alpha                       ( ln(1/T) - C3 )
    !
    !         r: reduced track length (l/l_0)
    !         l_0: initial, unannealed length
    !         t: Time of interest
    !         T: Temperature
    !         alpha: Empirical parameter fit to data
    !         beta: Empirical parameter fit to data
    !         C0: Empirical parameter fit to data
    !         C1: Empirical parameter fit to data
    !         C3: Empirical parameter fit to data
    !         C4: Empirical parameter fit to data
    !
    ! Inputs:
    !   nft_init:     number of fission tracks generated each fission event
    !   len_init:     nft_init-point double precision initial fission track length array (microns)
    !   nt:           number of time steps
    !   temp_celsius: n-point double precision temperature history array (Celsius)
    !   dep_km:       n-point double precision depth history array (km)
    !   dt_ma:        double precision timestep size (Ma)
    !   kinetic_par:  string indicating kinetic parameters (i.e., empirical parameters)
    !
    ! Outputs:
    !   ft_len:       (nft_init x nt) double precision fission track length array (microns)
    !----


    implicit none


    ! Arguments
    integer, intent(in) :: nft_init
    double precision, intent(in) :: len_init(nft_init)
    integer, intent(in) :: nt
    double precision, intent(in) :: temp_celsius(nt)
    double precision, intent(in) :: dep_km(nt)
    double precision, intent(in) :: dt_ma
    character(len=*), intent(in) :: kinetic_par
    double precision, intent(out) :: ft_len(nft_init,nt)


    ! Local variables
    integer :: i
    integer :: j
    double precision :: c0, c1, c2, c3, alpha, beta
    double precision :: dpar
    double precision :: r_mr_0
    double precision :: kappa
    double precision :: r(nt)
    double precision :: time
    double precision :: dt_seconds
    double precision :: temp_kelvin(nt)
    character(len=32) :: model_type
    logical :: keepTracks


    ! Select parameters for annealing
    model_type = 'none-specified'
    if (kinetic_par.eq.'green-et-al-1986') then
        c0 = -18.954d0
        c1 = 8.2385d-4
        c2 = -28.143d0
        c3 = 1.1217d-11
        alpha = -0.27951d0
        beta = -1.7910d0
        model_type = 'fanning-arrhenius'
    elseif (kinetic_par.eq.'ketcham-et-al-1999-all-fanning-arrhenius') then
        c0 = -11.053d0
        c1 = 3.8964d-4
        c2 = -17.842d0
        c3 = 6.7674d-4
        alpha = -0.14840d0
        beta = -8.7627d0
        model_type = 'fanning-arrhenius'
    elseif (kinetic_par.eq.'ketcham-et-al-1999-all-fanning-curvilinear') then
        c0 = -26.039d0
        c1 = 0.53168d0
        c2 = -62.319d0
        c3 = -7.8935d0
        alpha = -0.20196d0
        beta = -7.4224d0
        model_type = 'fanning-curvilinear'
        dpar = 2.50d0
    elseif (kinetic_par.eq.'ketcham-et-al-1999-all-fanning-curvilinear-c') then
        c0 = -19.844d0
        c1 = 0.38951d0
        c2 = -51.523d0
        c3 = -7.6423d0
        alpha = -0.12327d0
        beta = -11.988d0
        model_type = 'fanning-curvilinear'
        dpar = 1.75d0
    endif


    ! Convert units to Kelvin and seconds
    temp_kelvin = temp_celsius + 273.0d0
    dt_seconds = dt_ma * 1d6 * 365d0 * 24.0d0 * 60.0d0 * 60.0d0 ! SHOULD YEAR BE 365.25?


    ! Set constants for fission track parameters in different apatites
    r_mr_0 = 1.0d0 - exp(0.647d0*(dpar-1.75d0)-1.834d0)
    kappa = 1.0d0 - r_mr_0


    ! Initialize reduced track length array
    r = 0.0d0

    ! Calculate cumulative annealing backwards in time (Crowley, 1993)
    do i = nt,1,-1

        ! Calculate equivalent annealing time for current reduced length, temperature
        if (i.eq.nt) then
            time = exp(c3)
        else
            time = (((1.0d0-r(i+1)**beta)/beta)**alpha-1.0d0)/alpha
            if (model_type.eq.'fanning-arrhenius') then
                time = exp(((time-c0)/c1)*(1.0d0/temp_kelvin(i)-c3)+c2)
            elseif (model_type.eq.'fanning-curvilinear') then
                time = exp(((time-c0)/c1)*(log(1.0d0/temp_kelvin(i))-c3)+c2)
            else
                write(0,*) 'no model type called "',trim(model_type),'"'
                call error_exit(1)
            endif
        endif

        ! Add timestep duration to equivalent time
        time = time + dt_seconds

        ! Update reduced length
        if (model_type.eq.'fanning-arrhenius') then
            r(i) = c0 + c1*((log(time)-c2)/(1.0d0/temp_kelvin(i)-c3))
        elseif (model_type.eq.'fanning-curvilinear') then
            r(i) = c0 + c1*((log(time)-c2)/(log(1.0d0/temp_kelvin(i))-c3))
        else
            write(0,*) 'no model type called "',trim(model_type),'"'
            call error_exit(1)
        endif
        r(i) = (1.0d0 + alpha*r(i))**(1/alpha)
        r(i) = (1.0d0 - beta*r(i))**(1.0d0/beta)
    enddo


    ! Calculate the present-day lengths of fission tracks generated at each timestep
    do j = 1,nft_init
        do i = 1,nt

        ! Multiply initial track length by reduced length factor
        if (r(i)-r_mr_0.le.0.0d0) then
            ft_len(j,i) = 0.0d0
        else
            ft_len(j,i) = len_init(j) * ((r(i)-r_mr_0)/(1.0d0-r_mr_0))**(kappa)
        endif

        ! Set negative lengths to zero
            if (ft_len(j,i).lt.0.0d0) then
                ft_len(j,i) = 0.0d0
            endif

        enddo
    enddo


    ! Remove tracks from time period before mineral actually formed
    keepTracks = .false. ! At start of run, assume mineral has not formed and do not keep tracks
    do i = 1,nt
        if (.not.keepTracks) then
            if (dep_km(i).le.0.0d0) then
                keepTracks = .true.     ! Below surface, mineral exists, keep tracks from here on
            else
                ft_len(:,i) = 0.0d0     ! Above surface, remove tracks from this time
            endif
        endif
    enddo

    return
    end subroutine generate_fts_ketcham_et_al_1999


!--------------------------------------------------------------------------------------------------!



    subroutine segment_fts_carlson_1990( &
        nbins,                           &
        binwid,                          &
        len_min,                         &
        hist                             &
    )

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
    !     s:         empirically derived slope from fraction segmented vs. length plot (Fig. 8b)
    !
    ! Inputs:
    !   nbins:      number of histogram bins
    !   binwid:     histogram bin width (microns)
    !   len_min:    minimum length in histogram (microns)
    !
    ! Outputs:
    !   hist:       input track histogram with segmentation effects added (overwrites input array)
    !----


    implicit none


    ! Arguments
    integer, intent(in) :: nbins
    double precision, intent(in) :: binwid
    double precision, intent(in) :: len_min
    integer, intent(inout) :: hist(nbins)


    ! Local variables
    integer :: i
    integer :: j
    integer :: n
    integer :: ct
    integer :: ibin
    double precision :: len
    double precision, parameter :: s = 0.10d0 ! empirically determined (Carlson, 1990, Fig. 8b)
    double precision :: fraction_segmented
    double precision :: arg
    double precision :: x
    integer, allocatable :: seed(:)           ! random number variables
    integer :: epoch_time
    double precision :: ran



    ! Initialize a random number generator (for determining track segmentation points)
    if (.not.allocated(seed)) then
        epoch_time = time()
        call random_seed(SIZE=n)
        allocate(seed(n))
        seed(1:n) = epoch_time
        call random_seed(PUT=seed)
    endif



    ! For each length bin...
    do i = 1,nbins

        ! Fission track length in the bin
        len = dble(i-1)*binwid + len_min


        ! Calculate fraction of tracks that are segmented at this length (Carlson, 1990)
        if (len.lt.12.0d0) then
            fraction_segmented = s * (12.0d0 - len)
        else
            fraction_segmented = 0.0d0
            cycle
        endif


        ! For each track in this bin, determine whether it is segmented or not
        ! If the track is segmented, divide it into two random length tracks
        ct = 0
        do j = 1,hist(i)

            arg = dble(j)*fraction_segmented - dble(ct)   ! test for segmenting this track


            if (arg.ge.1.0d0) then

                !***** SEGMENT TRACK *****!

                call random_number(ran)                   !---- RANDOM SEGMENTATION SITE

                x = len*ran                               !---- SEGMENT 1
                ibin = int((x+0.5d0*binwid)/binwid)       ! which bin is segmented track 1 in?
                if (1.le.ibin.and.ibin.le.nbins) then
                    hist(ibin) = hist(ibin) + 1           ! add segmented track 1 to bin
                endif

                x = len - x                               !---- SEGMENT 2
                ibin = int((x+0.5d0*binwid)/binwid)       ! which bin is segmented track 2 in?
                if (1.le.ibin.and.ibin.le.nbins) then
                    hist(ibin) = hist(ibin) + 1           ! add segmented track 2 to bin
                endif

                hist(i) = hist(i) - 1                     !---- REMOVE INITIAL TRACK

                ! Update test count
                ct = ct + 1

            endif


        enddo


        ! Set negative counts to zero
        if (hist(i).lt.0) then
            hist(i) = 0
        endif


    enddo

    return

    end subroutine segment_fts_carlson_1990



!--------------------------------------------------------------------------------------------------!



    subroutine correct_fts_etching_userbias_willett_1997( &
        nbins,                                            &
        binwid,                                           &
        len_min,                                          &
        len0,                                             &
        hist                                              &
    )

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
    !   len0:       initial average fission track length (microns)
    !
    ! Outputs
    !   hist:       input track histogram with bias effects added (overwrites input array)
    !----


    implicit none


    ! Arguments
    integer, intent(in) :: nbins
    double precision, intent(in) :: binwid
    double precision, intent(in) :: len_min
    double precision, intent(in) :: len0
    integer, intent(inout) :: hist(nbins)


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


        ! Correct histogram for etching and user bias based on Willett (1997) Equation 4 fitting
        ! Green et al. (1988) density-versus-track length data
        if (len.le.a) then
            ! No short tracks
            hist(i) = 0
        elseif (len.le.b) then
            ! Reduce number of medium length tracks
            ratio = len/len0   ! ratio of bin length to average initial length
            hist(i) = int( dble(hist(i)) * ((2.862d0*ratio)-1.2104d0) )
        else
            ! Keep number of long tracks
            hist(i) = hist(i)
        endif


        ! Set negative counts to zero
        if (hist(i).lt.0) then
            hist(i) = 0
        endif

    enddo


    return

    end subroutine correct_fts_etching_userbias_willett_1997





!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!----------------------------------- FISSION TRACK AGE ROUTINES -----------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!



    subroutine calc_ft_retention_age( &
        nbins,                        &
        hist,                         &
        nft_init,                     &
        dt_ma,                        &
        age_ma                        &
    )

    !----
    ! Fission track "retention age" (Willett, 1997) is determined by counting the number of tracks.
    ! "It is obtained from the annealing model by noting what time tracks first form and survive to
    ! the present time." In other words, the more fission tracks, the longer the sample has been
    ! accumulating them at low enough temperatures to retain them.
    !
    !   # FTs = (Time Accumulating FTs) * (# FTs Produced Per Timestep) / (Timestep Duration)
    !
    !         ...rearranging...
    !
    !   age_ma = ntracks * dt_ma / nft_init
    !
    ! Inputs:
    !   nbins:      number of histogram bins
    !   hist:       histogram of track lengths
    !   nft_init:   number of fission tracks produced each timestep
    !   dt_ma:      time step (Ma)
    !
    ! Outputs
    !   age_ma:     retention age (Ma)
    !----


    implicit none


    ! Arguments
    integer, intent(in) :: nbins
    integer, intent(in) :: hist(nbins)
    integer, intent(in) :: nft_init
    double precision, intent(in) :: dt_ma
    double precision, intent(out) :: age_ma


    ! Local variables
    integer :: i
    integer :: ntracks



    ! Calculate total number of tracks in fission track length distribution
    ntracks = 0
    do i = 1,nbins
        ntracks = ntracks + hist(i)
    enddo


    ! Calculate retention age
    age_ma = ntracks * dt_ma / nft_init

    return

    end subroutine calc_ft_retention_age



!--------------------------------------------------------------------------------------------------!



    subroutine calc_ft_age( &
        nbins,              &
        binwid,             &
        len_min,            &
        hist,               &
        nft_init,           &
        avg_len_init,       &
        dt_ma,              &
        age_ma              &
    )

    !----
    ! The "fission track age" (Willett, 1997) is based on integrating the annealing process over
    ! time, weighting the tracks generated at each time component by their final mean length. In
    ! this age calculation, longer fission tracks (less annealed) contribute more to the age than
    ! shorter fission tracks (more annealed).
    !
    ! Willett (1997) Equation 7 describes the fission track age for continuously generated fission
    ! tracks in a model over time:
    !
    !            1    __
    !     age = --- * \  w_i * dt_i
    !           rst   /_
    !                  i
    !
    !         age: fission track age
    !         rst: factor for reducing age of standard in zeta dating method (0.893 or 1/0.893?)
    !         w_i: weighting factor for timestep i, p_i/p_0 ~ l_i/l_0
    !         dt_i: model timestep size
    !
    !         p_i: final track density from tracks produced at timestep i
    !         p_0: initial, unannealed track density
    !         l_i: final average track length from tracks produced at timestep i
    !         l_0: initial, unannealed track length
    !
    ! Assuming constant timestep size, dt, the above equation becomes:
    !
    !            1         __
    !     age = --- * dt * \  w_i
    !           rst        /_
    !                       i
    !
    ! The weighting factor is proportional to the reduced track density and therefore the track
    ! lengths. Expanding that relation to include all tracks generated at each timestep:
    !
    !            l_i       total_track_length_i
    !     w_i = ----- = ----------------------------
    !            l_0     total_initial_track_length
    !
    ! The age formula thus becomes:
    !
    !     age = 1/rst * dt * total_track_length / total_initial_track_length
    !
    ! Inputs:
    !   nbins:             number of histogram bins
    !   binwid:            histogram bin width (microns)
    !   len_min:           minimum fission track bin length (microns)
    !   hist:              histogram of track lengths
    !   nft_init:          number of fission tracks produced each timestep
    !   avg_len_init:      initial average track length (microns)
    !   dt_ma:             time step (Ma)
    !
    ! Outputs
    !   age_ma:            fission track age (Ma)
    !----


    implicit none


    ! Arguments
    integer, intent(in) :: nbins
    double precision, intent(in) :: binwid
    double precision, intent(in) :: len_min
    integer, intent(in) :: hist(nbins)
    integer, intent(in) :: nft_init
    double precision, intent(in) :: avg_len_init
    double precision, intent(in) :: dt_ma
    double precision, intent(out) :: age_ma


    ! Local variables
    double precision :: total_length_microns
    double precision, parameter :: rst = 0.893d0 ! Ratio between observed and induced FT length in
                                                 ! Durango apatite standard



    ! Calculate total number of tracks in fission track length distribution
    total_length_microns = total_len(nbins,binwid,len_min,hist)


    ! Calculate fission track age
    age_ma = total_length_microns * dt_ma / dble(nft_init) / avg_len_init !/ rst


    return

    end subroutine calc_ft_age




!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------- HISTOGRAM ROUTINES ---------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

    subroutine get_histogram_parameters(n, array, binwid, xmin, xmax, nbin)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: array(n)
    double precision, intent(in) :: binwid
    integer :: nbin
    double precision :: xmin
    double precision :: xmax
    xmin = dble(floor(minval(array)) - 1)
    xmax = dble(ceiling(maxval(array)) + 1)
    nbin = nint((xmax-xmin+1.0d0)/binwid)
    return
    end subroutine

!--------------------------------------------------------------------------------------------------!

    subroutine load_histogram(n, array, binwid, xmin, nbin, hist)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: array(n)
    double precision, intent(in) :: binwid
    integer, intent(in) :: nbin
    double precision, intent(in) :: xmin
    integer, intent(out) :: hist(nbin)
    integer :: i
    integer :: ibin
    hist = 0
    do i = 1,n
        ibin = nint((array(i) - xmin + binwid)/binwid) ! bin value is at center
        if (0.lt.ibin .and. ibin.le.nbin) then
            hist(ibin) = hist(ibin) + 1
        endif
    enddo
    return
    end subroutine

!--------------------------------------------------------------------------------------------------!

    function total_len(n,binwid,len_min,hist)
    implicit none
    double precision :: total_len
    integer :: n
    double precision :: binwid
    double precision :: len_min
    integer :: hist(n)
    integer :: i
    total_len = 0.0d0
    do i = 1,n
        total_len = total_len + dble(hist(i))*(len_min+binwid*dble(i-1))
    enddo
    return
    end function


end module fission_track