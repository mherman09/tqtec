cmake_minimum_required(VERSION 3.12)
project(TQTec VERSION 2022.11)
enable_language(C Fortran)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wextra -Wunused -Wconversion -fcheck=bounds,do,mem,pointer,recursion -fbacktrace")

add_library(tqteclib
            src/diffusion_module.f90
            src/error_exit.c
            src/fission_track_module.f90
            src/apatite_helium_module.f90
            src/radiogenic_helium_module.f90
)

add_executable(tqtec src/tqtec.f90)
add_dependencies(tqtec tqteclib)
target_link_libraries(tqtec tqteclib)

add_executable(minage src/minage.f90)
add_dependencies(minage tqteclib)
target_link_libraries(minage tqteclib)