# Need 3.7 for submodule support
cmake_minimum_required(VERSION 3.7)
project( proj Fortran )

#=============================================================================80
# CMake macros
#=============================================================================80
include( ${CMAKE_SOURCE_DIR}/cmake/macros/add_exe_sources.txt  )
include( ${CMAKE_SOURCE_DIR}/cmake/macros/add_lib_sources.txt  )

# set( _common_libraries "lib" )
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules")

#=============================================================================80
# Set output directories
#=============================================================================80
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/lib )

set(CMAKE_GCC_FLAGS "" )

# Custom debug build flags
set(
  GCC_USER_DEBUG_FLAGS
  "-g -Wall -Wextra -fimplicit-none -fcheck=all \
   -ffpe-trap=zero,underflow,overflow,invalid --coverage -fbacktrace \
   -fstack-protector-all -fstack-check -Wimplicit-procedure \
   -Wno-use-without-only -Wno-unused-dummy-argument"
)

# Custom release build flags
set( GCC_USER_RELEASE_FLAGS "-march=x86-64" )

# Standard cmake build flags
set( CMAKE_Fortran_LINK_FLAGS_DEBUG "" )
set( CMAKE_Fortran_LINK_FLAGS_RELEASE "" )

# Append custom build flags to cmake flags
set(
  CMAKE_Fortran_FLAGS_DEBUG
  "${CMAKE_Fortran_FLAGS_DEBUG} ${USER_DEBUG_FLAGS}"
)
set(
  CMAKE_Fortran_FLAGS_RELEASE
  "${CMAKE_Fortran_FLAGS_RELEASE} ${USER_RELEASE_FLAGS}"
)

#=============================================================================80
# Add source tree and included folders
#=============================================================================80
add_subdirectory(src)

#=============================================================================80
# Add library and executables
#=============================================================================80
add_library( lib ${LIB_SRCS} )
target_link_libraries( lib )
set_target_properties( lib PROPERTIES
                            COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}"
                            LINK_FLAGS    "${CMAKE_FORTRAN_LINK_FLAGS}" )

add_executable( proj ${proj_sources} )
target_link_libraries( proj lib )
set_target_properties( proj PROPERTIES
                            COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}"
                            LINK_FLAGS    "${CMAKE_FORTRAN_LINK_FLAGS}" )
