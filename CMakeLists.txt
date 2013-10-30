cmake_minimum_required( VERSION 2.8 )

project( FreeON Fortran C CXX )

set( PACKAGE_VERSION 1.0.7 )

add_subdirectory( Modules )
add_subdirectory( FreeON )
add_subdirectory( tools )

option( ENABLE_GIT_TAG
  "Adds the HEAD git tag to the binary, which might be useful when using
FreeON directly from git." )

set( CMAKE_Fortran_FLAGS "-j" )
try_compile( FC_LONG_LINE
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/cmake-tests/fortran_long_line.F90)

configure_file( ${CMAKE_SOURCE_DIR}/MondoConfig.h.in
  ${CMAKE_BINARY_DIR}/MondoConfig.h )
