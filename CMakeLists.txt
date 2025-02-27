cmake_minimum_required(VERSION 3.12)
project(TQTec VERSION 2024.10)
enable_language(C Fortran)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wextra -Wunused -Wconversion -fcheck=bounds,do,mem,pointer,recursion -fbacktrace")


# TQTec libraries
add_library(
    tqteclib
    src/diffusion_module.f90
    src/error_exit.c
    src/fission_track_module.f90
    src/physical_constants_module.f90
    src/radiogenic_helium_module.f90
)

if (NOT "${SUPERLU_DIR}" STREQUAL "" )
    find_package(BLAS)
    set(extra_libraries ${SUPERLU_DIR}/lib/libsuperlu_fortran.a ${SUPERLU_DIR}/lib/libsuperlu.a ${BLAS_LIBRARIES})
    add_compile_definitions(COMPILE_WITH_SUPERLU)
else()
    set(extra_libraries "")
    add_compile_definitions(COMPILE_WITHOUT_SUPERLU)
endif()


# TQTec executables
add_executable(tqtec src/tqtec.f90 src/tqtec_io.f90 src/tqtec_actions.f90)
add_dependencies(tqtec tqteclib)
target_link_libraries(tqtec tqteclib)
target_compile_definitions(tqtec PUBLIC COMPILE_TQTEC)
install(TARGETS tqtec DESTINATION bin)

add_executable(tqclim src/tqtec.f90 src/tqtec_io.f90 src/tqtec_actions.f90)
add_dependencies(tqclim tqteclib)
target_link_libraries(tqclim tqteclib)
target_compile_definitions(tqclim PUBLIC COMPILE_TQCLIM)
install(TARGETS tqclim DESTINATION bin)

add_executable(readtqtec src/readtqtec.f90)
add_dependencies(readtqtec tqteclib)
target_link_libraries(readtqtec tqteclib)
install(TARGETS readtqtec DESTINATION bin)

add_executable(minage src/minage.f90)
add_dependencies(minage tqteclib)
target_link_libraries(minage tqteclib ${extra_libraries})
install(TARGETS minage DESTINATION bin)

add_executable(petro src/petro.f90)
add_dependencies(petro tqteclib)
target_link_libraries(petro tqteclib)
install(TARGETS petro DESTINATION bin)

add_executable(build_temp_history src/build_temp_history.f90)
install(TARGETS build_temp_history DESTINATION bin)

# add_executable(ftage src/ftage.f90)
# install(TARGETS ftage DESTINATION bin)


# TQTec scripts
list(
    APPEND shell_scripts
    do_full_example.sh
    do_minage_example.sh
    do_tqtec_example.sh
    plot_minage_aft_final.sh
    plot_tqtec_all.sh
    plot_tqtec_closure.sh
    plot_tqtec_dep_vs_time.sh
    plot_tqtec_hf_vs_time.sh
    plot_tqtec_temp_contours.sh
    plot_tqtec_temp_vs_dep.sh
    plot_tqtec_temp_vs_time.sh
    plot_tqtec_tti_vs_time.sh
)
foreach(s ${shell_scripts})
    add_custom_target(${s}-script ALL DEPENDS ${s})
    add_custom_command(
        OUTPUT ${s}
        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/scripts/${s} ${CMAKE_CURRENT_BINARY_DIR}/${s}
        MAIN_DEPENDENCY ${CMAKE_SOURCE_DIR}/scripts/${s}
        DEPENDS tqtec readtqtec minage petro
        COMMENT "Building target ${s}"
    )
    install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${s} DESTINATION bin)
endforeach()
