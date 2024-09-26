!--------------------------------------------------------------------------------------------------!
! Petro                                                                                            !
!                                                                                                  !
! Authors:                                                                                         !
!     - Matt Herman (current Modern Fortran version, i.e., what you are looking at right now!)     !
!     - Matt Legg (PSU MS thesis)                                                                  !
!                                                                                                  !
! Petro calculates petroleum observables for a temperature-(depth-)time history. Currently, petro  !
! supports calculating:                                                                            !
!     - Lopatin's Time Temperature Index of Maturation (TTI)                                       !
!     - Vitrinite Reflectance                                                                      !
!     - Types I, II, and III Organic Keragen Production                                            !
!                                                                                                  !
! Petro is designed to work with TQTec and readTQTec output files, but will work with any          !
! temperature history in the format:                                                               !
!                                                                                                  !
!     total_time(Ma) timestep(Ma) nsteps                                                           !
!     temp_1                                                                                       !
!     temp_2                                                                                       !
!     temp_3                                                                                       !
!       :                                                                                          !
!     temp_nsteps                                                                                  !
!                                                                                                  !
! References                                                                                       !
! Lopatin, N.V. (1971). Temperature and geologic time as factors in coalification (in Russian).    !
!     Akad. Nauk SSSR Izv. Ser. Geol., 3, 95-106.                                                  !
! Waples, D.W. (1980). Time and temperature in petroleum formation: Application of Lopatin's       !
!     method to petroleum exploration. AAPG Bulletin, 64(6), 916-926.                              !
!--------------------------------------------------------------------------------------------------!

module petro

    ! Input files
    character(len=512) :: readtqtec_temp_file
    character(len=512) :: readtqtec_dep_file

    ! Output files
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

    ! TTI variables
    character(len=512) :: tti_file

    ! Vitrinite reflectance variables
    character(len=512) :: vitrinite_file

end module petro

!==================================================================================================!

program main

use petro, only:         &
    readtqtec_temp_file, &
    readtqtec_dep_file,  &
    isOutputDefined,     &
    tti_file,            &
    vitrinite_file


implicit none




! Parse command line arguments
call gcmdln()




! Input/output error checks

! Is input temperature history defined?
if (readtqtec_temp_file.eq.'') then
    call usage('petro: no readtqtec output temperature file defined')
endif

! Is input depth history defined?
if (readtqtec_dep_file.eq.'') then
    write(0,*) 'petro: no readtqtec output depth file defined - setting all depths to 0 km'
endif

! Is an output defined?
if (.not.isOutputDefined) then
    call usage('petro: no output defined')
endif




! Read temperature history
call read_temp_history()



! Read depth history
call read_dep_history()




!***************************************************!
!********** CALCULATE PETROLEUM VARIABLES **********!
!***************************************************!

! Time-temperature index of maturation
if (tti_file.ne.'') then
    call calc_tti()
endif

! Vitrinite reflectance
if (vitrinite_file.ne.'') then
    call calc_vitrinite_reflectance()
endif



end program main




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

use petro, only:         &
    readtqtec_temp_file, &
    nhorizons,           &
    nhorizons_max,       &
    time_ma,             &
    dt_ma,               &
    ntimes,              &
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
    write(0,*) 'petro: error reading first line of tqtec output temperature file'
    call error_exit(1)
endif


! Extract number of horizons tracked in the thermal model from first line of temperatures
nhorizons = 1
read(21,'(A)',iostat=ios) line
if (ios.ne.0) then
    write(0,*) 'petro: error reading second line of tqtec output temperature file'
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
        write(0,*) 'petro: error reading line',i,' of readtqtec output temperature file'
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

use petro, only:        &
    readtqtec_dep_file, &
    nhorizons,          &
    nhorizons_max,      &
    time_ma,            &
    dt_ma,              &
    ntimes,             &
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
    write(0,*) 'petro: error reading first line of tqtec output depth file'
    call error_exit(1)
endif
if (ntimes_dep.ne.ntimes) then
    call usage('petro: number of timesteps differs between temp and dep files')
endif


! Extract number of horizons tracked in the thermal model from first line of temperatures
nhorizons_dep = 1
read(20,'(A)',iostat=ios) line
if (ios.ne.0) then
    write(0,*) 'petro: error reading second line of tqtec output depth file'
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
    call usage('petro: number of tracked horizons differs between temp and dep files')
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
        write(0,*) 'petro: error reading line',i,' of readtqtec output depth file'
        call error_exit(1)
    endif
enddo


! Close the depth file
close(20)


return
end subroutine





!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!-------------------------- PETROLEUM CALCULATION SUBROUTINES -------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!


subroutine calc_tti()

use petro, only: &
    nhorizons,          &
    dt_ma,              &
    ntimes,             &
    temp_celsius_array, &
    dep_km_array,       &
    tti_file

implicit none


! Local variables
integer :: i
integer :: j
double precision, allocatable :: tti(:)
double precision, allocatable :: tti_exp(:)
character(len=32) :: fmt_string


write(*,*) 'petro: calculating time-temperature index of maturity (TTI)'


! Allocate memory to TTI arrays
if (allocated(tti)) then
    deallocate(tti)
endif
allocate(tti(nhorizons))
if (allocated(tti_exp)) then
    deallocate(tti_exp)
endif
allocate(tti_exp(nhorizons))


! Open output file
open(unit=31,file=tti_file,status='unknown')


! Header: time_before_present(Ma)  sample_rate(Ma)  number_of_points
write(31,'(2F10.3,I6)') -ntimes*dt_ma, dt_ma, ntimes

! Calculate and print TTI for each horizon at each timestep
write(fmt_string,'("("I6,"E14.6)")') nhorizons
do i = 1,ntimes
    tti_exp = 0.1d0*(temp_celsius_array(:,i)-105.0d0)
    if (i.eq.1) then
        tti = dt_ma * 2**tti_exp
    else
        tti = tti + dt_ma * 2**tti_exp
    endif
    write(31,fmt_string) (tti(j),j=1,nhorizons)
enddo

close(31)

return
end subroutine


!--------------------------------------------------------------------------------------------------!


subroutine calc_vitrinite_reflectance()

use petro, only:        &
    nhorizons,          &
    dt_ma,              &
    ntimes,             &
    temp_celsius_array, &
    dep_km_array,       &
    vitrinite_file

implicit none


! Local variables
integer :: i
integer :: j
double precision, allocatable :: tti(:)
double precision, allocatable :: tti_exp(:)
double precision, allocatable :: vr0(:)
character(len=32) :: fmt_string


write(*,*) 'petro: calculating vitrinite reflectance'


! Allocate memory to TTI arrays
if (allocated(tti)) then
    deallocate(tti)
endif
allocate(tti(nhorizons))
if (allocated(tti_exp)) then
    deallocate(tti_exp)
endif
allocate(tti_exp(nhorizons))
if (allocated(vr0)) then
    deallocate(vr0)
endif
allocate(vr0(nhorizons))



! Open output file
open(unit=32,file=vitrinite_file,status='unknown')


! Header: time_before_present(Ma)  sample_rate(Ma)  number_of_points
write(32,'(2F10.3,I6)') -ntimes*dt_ma, dt_ma, ntimes

! Calculate and print TTI for each horizon at each timestep
write(fmt_string,'("("I6,"E14.6)")') nhorizons
do i = 1,ntimes
    tti_exp = 0.1d0*(temp_celsius_array(:,i)-105.0d0)
    if (i.eq.1) then
        tti = dt_ma * 2**tti_exp
    else
        tti = tti + dt_ma * 2**tti_exp
    endif
    vr0 = (tti/41.0d0) ** 0.22d0
    write(31,fmt_string) (vr0(j),j=1,nhorizons)
enddo

close(32)


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

use petro, only:         &
    readtqtec_temp_file, &
    readtqtec_dep_file,  &
    isOutputDefined,     &
    tti_file,            &
    vitrinite_file

implicit none

integer :: i
integer :: narg
integer :: ios
character(len=512) :: tag
logical :: fileExists


! Initialize variables
readtqtec_temp_file = ''
readtqtec_dep_file = ''
isOutputDefined = .false.
tti_file = ''
vitrinite_file = ''


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
    elseif (trim(tag).eq.'-tti') then
        i = i + 1
        call get_command_argument(i,tti_file,status=ios)
        isOutputDefined = .true.

    elseif (trim(tag).eq.'-vitrinite') then
        i = i + 1
        call get_command_argument(i,vitrinite_file,status=ios)
        isOutputDefined = .true.


    else
        call usage('petro: no option '//trim(tag))
    endif

    if (ios.ne.0) then
        call usage('petro: error parsing "'//trim(tag)//'" flag arguments')
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
write(0,*) 'Usage: petro -temp READTQTEC_TEMP_FILE [-dep READTQTEC_DEP_FILE] [...options...]'
write(0,*)
write(0,*) '-temp READTQTEC_TEMP_FILE  Temperature history output from readtqtec'
write(0,*) '-dep READTQTEC_DEP_FILE    Depth history (default: track age from beginning of model)'
write(0,*) '-tti TTI_FILE              Lopatin time-temperature index of maturation'
write(0,*) '-vitrinite VIT_FILE        Vitrinite reflectance'
write(0,*) '-advanced                  See advanced options'
write(0,*)
call error_exit(1)
stop
end subroutine
