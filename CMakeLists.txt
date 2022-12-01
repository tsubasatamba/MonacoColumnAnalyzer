####### CMakeLists.txt
cmake_minimum_required(VERSION 3.0)

### Initial definition of cmake variables
set(CMAKE_INSTALL_PREFIX $ENV{HOME} CACHE PATH "install prefix")
set(CMAKE_BUILD_TYPE Release CACHE STRING "build type")
set(CMAKE_CXX_FLAGS_DEBUG "-g -W -Wall" CACHE STRING "CXX_FLAGS for debug")
set(CMAKE_C_FLAGS_DEBUG "-g -W -Wall" CACHE STRING "C_FLAGS for debug")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -W -Wall" CACHE STRING "CXX_FLAGS for release")
set(CMAKE_C_FLAGS_RELEASE "-O3 -W -Wall" CACHE STRING "C_FLAGS for release")
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

### Definition of project
project(BinaryAnalyzer)
add_definitions("-std=c++17")

### BOOST ###
find_package(Boost 1.56.0 REQUIRED COMPONENTS system filesystem)
set(BOOST_INC_DIR ${Boost_INCLUDE_DIRS})
set(BOOST_LIB_DIR ${Boost_LIBRARY_DIRS})
set(BOOST_LIB ${Boost_LIBRARIES})
message("-- BOOST_INC_DIR: ${BOOST_INC_DIR}")
message("-- BOOST_LIB_DIR: ${BOOST_LIB_DIR}")
message("-- BOOST_LIB: ${BOOST_LIB}")

### ANL ###
if(NOT DEFINED ANLNEXT_INSTALL)
  if(DEFINED ENV{ANLNEXT_INSTALL})
    set(ANLNEXT_INSTALL $ENV{ANLNEXT_INSTALL})
  else()
    set(ANLNEXT_INSTALL $ENV{HOME})
  endif()
endif(NOT DEFINED ANLNEXT_INSTALL)
set(ANLNEXT_INC_DIR ${ANLNEXT_INSTALL}/include)
set(ANLNEXT_LIB_DIR ${ANLNEXT_INSTALL}/lib)
set(ANLNEXT_LIB ANLNext)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${ANLNEXT_LIB_DIR}/anlnext)
message("-- ANLNEXT_INSTALL = ${ANLNEXT_INSTALL}")

### ROOT ###
set(ROOTSYS $ENV{ROOTSYS})
list(APPEND CMAKE_PREFIX_PATH ${ROOTSYS})
find_package(ROOT REQUIRED)
set(ROOT_INC_DIR ${ROOT_INCLUDE_DIRS})
set(ROOT_LIB_DIR ${ROOT_LIBRARY_DIR})
set(ROOT_LIB ${ROOT_LIBRARIES})
list(APPEND CMAKE_INSTALL_RPATH ${ROOT_LIBRARY_DIR})
message("-- ROOTSYS = ${ROOTSYS}")
message("-- ROOT_INC_DIR = ${ROOT_INC_DIR}")
message("-- ROOT_LIB_DIR = ${ROOT_LIB_DIR}")
message("-- ROOT libraries = ${ROOT_LIB}")
add_definitions(-DUSE_ROOT)
add_definitions(-DCS_BASIC2)

### Geant4 ###
find_package(Geant4 REQUIRED)
set(G4_INC_DIR ${Geant4_INCLUDE_DIRS})
set(G4_LIB ${Geant4_LIBRARIES})
message("-- Geant4 libraries: ${G4_LIB}")

### ANLGeant4
set(ANLG4_INC_DIR $ENV{HOME}/include/anlgeant4)
set(ANLG4_LIB_DIR $ENV{HOME}/lib)
set(ANLG4_LIB ANLGeant4)

### FITSIO ###
find_path(CFITSIO_INC_DIR
  NAMES fitsio.h
  PATHS /usr/local/include $ENV{HEADAS}/include)
find_library(CFITSIO_LIB
  NAMES cfitsio
  PATHS /usr/local/lib $ENV{HEADAS}/lib)

add_subdirectory(source)

