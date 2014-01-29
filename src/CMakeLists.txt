cmake_minimum_required (VERSION 2.6)

#SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
#SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_INSTALL_PREFIX})

# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
enable_language (Fortran)

add_executable(brktos brktos.f90 sjp_util.f90)
add_executable(rfacXN rfacXN2.f90 sjp_util.f90)
add_executable(aps aps.f90)
add_executable(frcc frcc.f90)
add_executable(hypro hypro.f90)
add_executable(sct sct.f90)
add_executable(sluv sluv.f90)

SET (EXES brktos rfacXN aps frcc hypro sct sluv)

install(TARGETS ${EXES} DESTINATION . COMPONENT binaries)
