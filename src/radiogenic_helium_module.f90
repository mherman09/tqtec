module radiogenic_helium

! Atomic masses of helium-producing radioactive elements (IUGS 1977)
double precision, parameter :: atomic_mass_th232 = 232.0380553d0
double precision, parameter :: atomic_mass_u235 = 235.043928d0
double precision, parameter :: atomic_mass_u238 = 238.0507826d0

! Decay constants of helium-producing radioactive elements (IUGS 1977)
double precision, parameter :: decay_th232 = 4.94750d-5/1.0d6/365.0d0/24.0d0/60.0d0/60.0d0  ! [1/s]
double precision, parameter :: decay_u235 =  9.84850d-4/1.0d6/365.0d0/24.0d0/60.0d0/60.0d0  ! [1/s]
double precision, parameter :: decay_u238 =  1.55125d-4/1.0d6/365.0d0/24.0d0/60.0d0/60.0d0  ! [1/s]



!==================================================================================================!
contains
!==================================================================================================!




subroutine calc_he_production_rate(mol_th232,mol_u235,mol_u238,volume,p)
!----
! Determine the radiogenic helium production rate [mol/m^3/s] in a volume, given the amount of
! radioactive thorium-232, uranium-235, and uranium-238 [moles]. The distribution of parent
! nuclides within the volume is assumed to be homogeneous.
!
! Inputs:
!   mol_th232:  moles of thorium-232
!   mol_u235:   moles of uranium-235
!   mol_u238:   moles of uranium-238
!   volume:     volume [m^3]
!
! Outputs:
!   p:          helium production rate [mol/m^3/s]
!----

implicit none

! Arguments
double precision :: mol_th232
double precision :: mol_u235
double precision :: mol_u238
double precision :: volume
double precision :: p

! Local variables
double precision :: mol_conc_th232
double precision :: mol_conc_u235
double precision :: mol_conc_u238


! Molar concentration (mol/m^3)
mol_conc_th232 = mol_th232/volume
mol_conc_u235 = mol_u235/volume
mol_conc_u238 = mol_u238/volume

! Helium production rate (mol/m^3/s)
! thorium-232 -> lead-208: 6 alpha particles (helium nuclei)
! uranium-235 -> lead-207: 7 alpha particles
! uranium-238 -> lead-206: 8 alpha particles
p = 8d0*mol_conc_u238*decay_u238 + 7d0*mol_conc_u235*decay_u235 + 6d0*mol_conc_th232*decay_th232

return
end subroutine




!--------------------------------------------------------------------------------------------------!




subroutine calc_u_th_he_age(mol_th232,mol_u235,mol_u238,mol_he4,age)
!----
! Determine the (U-Th)/He age based on the amount of parent uranium and thorium and daughter helium
! present in a sample.
!
! Inputs:
!   mol_th232:  moles of thorium-232
!   mol_u235:   moles of uranium-235
!   mol_u238:   moles of uranium-238
!   mol_he4:    moles of helium-4
!
! Outputs:
!   age:        Age [Ma]
!----

implicit none

! Arguments
double precision :: mol_th232
double precision :: mol_u235
double precision :: mol_u238
double precision :: mol_he4
double precision :: age

! Local variables
double precision :: t_ma
double precision :: dt_ma
double precision :: tmax_ma
double precision :: t_seconds
double precision :: dt_seconds
double precision :: tmax_seconds
double precision, parameter :: ma2s = 1.0d6*365.0d0*24.0d0*60.0d0*60.0d0
double precision :: f


! Initialize age search parameters
t_ma = 0.0d0              ! Initial age [Ma]
tmax_ma = 4.5d3           ! Maximum age [Ma]
dt_ma = 1.0d0             ! Starting age step size [Ma]

! Convert to seconds
t_seconds = t_ma*ma2s
tmax_seconds = tmax_ma*ma2s
dt_seconds = dt_ma*ma2s

! Solve for age by searching for zeros of helium production function
do while (t_seconds.le.tmax_seconds)

    ! Age function, f, is monotonically increasing, and its zero is the age
    f = 8.0d0*mol_u238*(exp(decay_u238*t_seconds)-1.0d0) + &
        7.0d0*mol_u235*(exp(decay_u235*t_seconds)-1.0d0) + &
        6.0d0*mol_th232*(exp(decay_th232*t_seconds)-1.0d0) - &
        mol_he4

    ! Crossed zero?
    if (f.gt.0.0d0) then
        t_seconds = t_seconds - dt_seconds
        dt_seconds = dt_seconds/10.0d0
        if (dt_seconds.lt.1d-6*ma2s) then
            exit
        endif
    endif

    t_seconds = t_seconds + dt_seconds

enddo

! Calculate age
age = t_seconds/ma2s

return
end subroutine


end module