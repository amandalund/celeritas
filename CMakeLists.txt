#---------------------------------*-CMake-*----------------------------------#
# Copyright 2020 UT-Battelle, LLC and other Celeritas Developers.
# See the top-level COPYRIGHT file for details.
# SPDX-License-Identifier: (Apache-2.0 OR MIT)
#----------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.10)
project(Celeritas VERSION 0.0.1 LANGUAGES CXX)
cmake_policy(VERSION 3.10...3.16)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
include(GNUInstallDirs)
include(CeleritasUtils)

#---------------------------------------------------------------------------##
# OPTIONS
#---------------------------------------------------------------------------##

# Tests
option(CELERITAS_BUILD_TESTS "Build Celeritas unit tests" ON)

# TPLs
option(CELERITAS_USE_CUDA "Enable GPU transport" ON)
option(CELERITAS_USE_MPI "Enable distributed memory parallelism" ON)
option(CELERITAS_USE_ROOT "Enable ROOT I/O" OFF)
option(CELERITAS_USE_VECGEOM "Enable VecGeom geometry" ON)

# Library
option(BUILD_SHARED_LIBS "Build shared libraries" ON)
if(NOT DEFINED CMAKE_INSTALL_RPATH_USE_LINK_PATH)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON CACHE STRING
      "Inform installed binaries of external library rpaths")
endif()
if(BUILD_SHARED_LIBS)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}" CACHE STRING
      "Inform installed binaries of internal library rpaths")
endif()
if(APPLE)
  option(CMAKE_MACOSX_RPATH "Support @rpath in install targets" ON)
endif()

# Build flags
option(CELERITAS_DEBUG "Enable runtime assertions" ON)
set(CMAKE_CUDA_FLAGS_DEBUG "-g -G" CACHE STRING "")
if(NOT DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Build type")
endif()

foreach(_LANG CUDA CXX)
  set(CMAKE_${_LANG}_STANDARD 14)
  set(CMAKE_${_LANG}_EXTENSIONS OFF)
endforeach()

#---------------------------------------------------------------------------##
# DEPENDENCIES
#---------------------------------------------------------------------------##

if(CELERITAS_USE_CUDA)
  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD 14)
endif()

if(CELERITAS_USE_MPI)
  find_package(MPI REQUIRED)
endif()

if(CELERITAS_USE_ROOT)
  find_package(ROOT REQUIRED)
endif()

if(CELERITAS_USE_VECGEOM)
  find_package(VecGeom REQUIRED)

  if(CELERITAS_USE_CUDA AND VECGEOM_CUDA_ARCH STREQUAL "none")
    message(FATAL_ERROR "The VecGeom installation at ${VECGEOM_INSTALL_DIR} "
      "was configured without CUDA, so CELERITAS_USE_CUDA cannot be used")
  endif()
endif()

add_subdirectory(external)

#---------------------------------------------------------------------------##
# LIBRARY
#---------------------------------------------------------------------------##

add_subdirectory(src)

#---------------------------------------------------------------------------##
# TESTS
#---------------------------------------------------------------------------##

if(CELERITAS_BUILD_TESTS)
  include(CTest)
  add_subdirectory(test)
endif()

#---------------------------------------------------------------------------##
# DEMO APP
#---------------------------------------------------------------------------##

# add_subdirectory(app)

#---------------------------------------------------------------------------##
