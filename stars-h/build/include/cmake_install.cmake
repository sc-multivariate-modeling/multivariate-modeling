# Install script for directory: /home/abdullsm/develop/stars-h/include

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/abdullsm/develop/stars-h/build/install_dir")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/home/abdullsm/develop/stars-h/include/starsh.h"
    "/home/abdullsm/develop/stars-h/include/starsh-mpi.h"
    "/home/abdullsm/develop/stars-h/include/starsh-starpu.h"
    "/home/abdullsm/develop/stars-h/include/starsh-mpi-starpu.h"
    "/home/abdullsm/develop/stars-h/include/starsh-spatial.h"
    "/home/abdullsm/develop/stars-h/include/starsh-spatial-gsl.h"
    "/home/abdullsm/develop/stars-h/include/starsh-electrostatics.h"
    "/home/abdullsm/develop/stars-h/include/starsh-electrodynamics.h"
    "/home/abdullsm/develop/stars-h/include/starsh-randtlr.h"
    "/home/abdullsm/develop/stars-h/include/starsh-minimal.h"
    "/home/abdullsm/develop/stars-h/include/starsh-constants.h"
    "/home/abdullsm/develop/stars-h/include/starsh-particles.h"
    "/home/abdullsm/develop/stars-h/include/starsh-cauchy.h"
    )
endif()

