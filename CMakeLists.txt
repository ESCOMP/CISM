# CMake file for CISM

cmake_minimum_required(VERSION 2.8.4)

OPTION (CISM_BUILD_CISM_DRIVER "Toggle to build cism_driver, on by default" ON)

OPTION (CISM_USE_GPTL_INSTRUMENTATION "Toggle to use GPTL instrumentation, on by default " ON)
OPTION (CISM_COUPLED "Toggle to build CISM for use with CESM, off by default" OFF)
# OPTION (CISM_USE_CISM_FRONT_END "Toggle to use cism_driver or cism_cesm_interface with cism_front_end, off by default" OFF)

OPTION (CISM_STATIC_LINKING "Toggle to set static linking for executables, off by default" OFF)
OPTION (CISM_FORCE_FORTRAN_LINKER "Toggle to force using a fortran linker for building executables, off by default" OFF)
OPTION (CISM_INCLUDE_IMPLICIT_LINK_LIBRARIES "Toggle to explicitly include the CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES on the link line, on by default" ON)

IF (CISM_BUILD_CISM_DRIVER)
  SET(CISM_USE_CISM_FRONT_END "ON")
ENDIF()

# set (or override) other options, if CISM_COUPLED is ON:
IF (CISM_COUPLED)
  SET(CISM_NO_EXECUTABLE "ON")
  SET(CISM_USE_CISM_FRONT_END "OFF")
  ADD_DEFINITIONS(-DCCSMCOUPLED)
ENDIF()


IF (NOT DEFINED CISM_BINARY_DIR)
   SET(CISM_BINARY_DIR ${CMAKE_BINARY_DIR})
ENDIF()

MESSAGE("Building in: ${CISM_BINARY_DIR}")

ENABLE_LANGUAGE(Fortran) 

IF (CISM_INCLUDE_IMPLICIT_LINK_LIBRARIES)
  # WJS (6-3-14) Until now, these Fortran_IMPLICIT_LINK_LIBRARIES were always
  # included. However, explicitly appending these implicit link libraries breaks
  # the build on yellowstone, for some reason. I suspect these implicit link
  # libraries are needed when you're linking with a C++ linker, which is not the
  # case for yellowstone-intel. In any case, I'm providing an option to exclude
  # these from the build, to allow the yellowstone-intel build to work.
  LIST(APPEND CISM_IMPLICIT_LINK_LIBS ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})

  #message("")
  #message("  CMake detected the following libraries for linking Fortran with C++ compiler:")
  #message("     ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES} ")
ENDIF()

MESSAGE(">> CISM_NETCDF_DIR set to : ${CISM_NETCDF_DIR}")

IF (NOT DEFINED CISM_NETCDF_DIR)
  MESSAGE(FATAL_ERROR "\nCISM Error: cmake must define CISM_NETCDF_DIR:
  (-D CISM_NETCDF_DIR=<netcdf_install_prefix>)!")
ENDIF()

FIND_PATH(CISM_NETCDFF_FOUND libnetcdff.a ${CISM_NETCDF_DIR}/lib)
#MESSAGE(" CISM_NETCDFF_FOUND = ${CISM_NETCDFF_FOUND}")

IF (${CISM_NETCDFF_FOUND} STREQUAL ${CISM_NETCDF_DIR}/lib )
   SET(CISM_NETCDF_LIBS "netcdff;netcdf" CACHE STRING "Netcdf Library Names(s)")
ELSE()
   SET(CISM_NETCDF_LIBS "netcdf" CACHE STRING "Netcdf Library Names(s)") 
ENDIF()
MESSAGE(">> CISM_NETCDF_LIBS Library(s) set to : ${CISM_NETCDF_LIBS}")

IF (DEFINED CISM_MPI_BASE_DIR)
  IF (NOT DEFINED CISM_MPI_LIB_DIR)
    SET(CISM_MPI_LIB_DIR ${CISM_MPI_BASE_DIR}/lib)
  ENDIF()
  IF (NOT DEFINED CISM_MPI_INC_DIR)
    SET(CISM_MPI_INC_DIR ${CISM_MPI_BASE_DIR}/include)
  ENDIF()
ENDIF()
MESSAGE(">> CISM_MPI_LIB_DIR set to : ${CISM_MPI_LIB_DIR}")
MESSAGE(">> CISM_MPI_INC_DIR set to : ${CISM_MPI_INC_DIR}")


IF (NOT DEFINED CMAKE_Fortran_MODULE_DIRECTORY)
  SET(CMAKE_Fortran_MODULE_DIRECTORY ${CISM_BINARY_DIR}/fortran_mod_files)
ENDIF()
INCLUDE_DIRECTORIES(${CMAKE_Fortran_MODULE_DIRECTORY})


PROJECT(CISM Fortran C)



# Copy a few needed files to the build directory:
MESSAGE(">> Calling utils/build/autocopy-to-build-dir")
EXECUTE_PROCESS(COMMAND  ${CISM_SOURCE_DIR}/utils/build/autocopy-to-build-dir
                ${CISM_SOURCE_DIR}
                WORKING_DIRECTORY ${CISM_BINARY_DIR}
                OUTPUT_FILE ${CISM_BINARY_DIR}/autocopy.log)

MESSAGE(">> see ${CISM_BINARY_DIR}/autocopy.log")
INCLUDE_DIRECTORIES(${CISM_BINARY_DIR}/fortran_autocopy_includes)

# End of setup and error checking

## Use CMAKE_Fortran_FLAGS to set flags; use CISM_Fortran_FLAGS to append flags
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${CISM_Fortran_FLAGS}")

ADD_DEFINITIONS(-DCISM_MPI)
ADD_DEFINITIONS(-DGLIMMER_MPI)
include_directories(${CISM_NETCDF_DIR}/include ${CISM_MPI_INC_DIR})

link_directories(${CISM_NETCDF_DIR}/lib ${CISM_MPI_LIB_DIR} ${CISM_SCI_LIB_DIR} 
                 ${CISM_HDF5_LIB_DIR} ) 

IF (CISM_USE_GPTL_INSTRUMENTATION)
 OPTION(CISM_GPTL_DIR "Optional path to installed gptl library")
 IF (CISM_GPTL_DIR)
   message(">> GPTL Library reqested, installed in ${CISM_GPTL_DIR}")
   ADD_DEFINITIONS(-DCESMTIMERS)
   include_directories(${CISM_GPTL_DIR})
   link_directories(${CISM_GPTL_DIR})
   SET(CISM_GPTL_LIB gptl)
 ELSE()
   message(">>GPTL Library not requested: can set CISM_GPTL_DIR to enable")
 ENDIF()
ENDIF()

#####Build All Fortran Sources #####

#IF (CISM_USE_CISM_FRONT_END)
   file(GLOB GLINT_SOURCES libglint/*.F90)
   file(GLOB GLAD_SOURCES libglad/*.F90)
   file(GLOB GLISSADE_SOURCES libglissade/*.F90)
#ENDIF()

file(GLOB FORTRANSOURCES
     libglimmer-solve/SLAP/*.f
     libglimmer-solve/*.F90
     libglimmer/*.F90 libglimmer/writestats.c
     libglide/*.F90
     cism_driver/eismint_forcing.F90
     cism_driver/testsfg.F90
     ${GLINT_SOURCES}
     ${GLAD_SOURCES}
     ${GLISSADE_SOURCES})

set (FORTRAN_AUTOGEN_SOURCES ${CISM_BINARY_DIR}/fortran_autogen_srcs/glide_io.F90
			     ${CISM_BINARY_DIR}/fortran_autogen_srcs/glint_io.F90
			     ${CISM_BINARY_DIR}/fortran_autogen_srcs/glint_mbal_io.F90
			     ${CISM_BINARY_DIR}/fortran_autogen_srcs/glad_io.F90
			     ${CISM_BINARY_DIR}/fortran_autogen_srcs/glad_mbal_io.F90
			     ${CISM_BINARY_DIR}/fortran_autogen_srcs/glimmer_vers.F90) 
message("Autogenerated CISM sources: ${FORTRAN_AUTOGEN_SOURCES}")

LIST(REMOVE_ITEM FORTRANSOURCES
   ${CISM_SOURCE_DIR}/libglimmer-solve/SLAP/dlapqc.f )



# Auto-generate the *_io.F90 files in the build directory:
MESSAGE(">> Calling utils/build/autogenerate-in-build-dir")

ADD_CUSTOM_COMMAND(
		OUTPUT fortran_autogen_srcs/glide_io.F90
		COMMAND ${CISM_SOURCE_DIR}/utils/build/autogenerate-in-build-dir
		${CISM_SOURCE_DIR} 
		DEPENDS ${CISM_SOURCE_DIR}/libglide/glide_vars.def 
		        ${CISM_SOURCE_DIR}/libglimmer/ncdf_template.F90.in
			${CISM_SOURCE_DIR}/utils/build/generate_ncvars.py
			${CISM_SOURCE_DIR}/utils/build/autogenerate-in-build-dir
)


# Auto-generate the *_io.F90 files in the build directory:
MESSAGE(">> Calling utils/build/autogen-for-glint-and-glad-in-build-dir")

ADD_CUSTOM_COMMAND(
		OUTPUT fortran_autogen_srcs/glint_io.F90
		       fortran_autogen_srcs/glint_mbal_io.F90
		       fortran_autogen_srcs/glad_io.F90
		       fortran_autogen_srcs/glad_mbal_io.F90
		COMMAND ${CISM_SOURCE_DIR}/utils/build/autogen-for-glint-and-glad-in-build-dir
                ${CISM_SOURCE_DIR}
		DEPENDS ${CISM_SOURCE_DIR}/libglint/glint_vars.def
			${CISM_SOURCE_DIR}/libglint/glint_mbal_vars.def
			${CISM_SOURCE_DIR}/libglad/glad_vars.def
			${CISM_SOURCE_DIR}/libglad/glad_mbal_vars.def
			${CISM_SOURCE_DIR}/libglimmer/ncdf_template.F90.in
			${CISM_SOURCE_DIR}/utils/build/generate_ncvars.py
			${CISM_SOURCE_DIR}/utils/build/autogen-for-glint-and-glad-in-build-dir
)

# Get autogenerated source files, and add them to the FORTRANSOURCES list:
# changed to always do this:
#IF (CISM_COUPLED)

# Remove old versions of autogenerated F90 files that may be sitting
# around in the source tree from an old cmake-based build or an
# autotools-based build. (Now the cmake-based build doesn't place
# these in the source tree.)
# Unlike SOURCEMODFILES, we just hard-code the files to remove for
# simplicity.
# This can be removed once we switch to consistently using this new
# cmake-based build.
  LIST(REMOVE_ITEM FORTRANSOURCES
     ${CISM_SOURCE_DIR}/libglide/glide_io.F90
     ${CISM_SOURCE_DIR}/libglide/glide_lithot_io.F90
     ${CISM_SOURCE_DIR}/libglint/glint_io.F90
     ${CISM_SOURCE_DIR}/libglint/glint_mbal_io.F90
     ${CISM_SOURCE_DIR}/libglimmer/glimmer_vers.F90 )
#ENDIF()


### (For CESM) Remove source files with names already in CISM_SOURCEMOD_DIR
OPTION(CISM_SOURCEMOD_DIR 
      "Path to SourceMod directory of F90 files to replace Glimmer files")

#MESSAGE("Fortran Source Files: ${FORTRANSOURCES}")

# Note that the following glob does NOT contain .cpp files, because
# those are built in a separate library - so for now, you can't put
# .cpp files in your sourceMod directory.
FILE(GLOB SOURCEMODFILES 
     ${CISM_SOURCEMOD_DIR}/*.F90
     ${CISM_SOURCEMOD_DIR}/*.F
     ${CISM_SOURCEMOD_DIR}/*.f90
     ${CISM_SOURCEMOD_DIR}/*.f
     ${CISM_SOURCEMOD_DIR}/*.c)

# MESSAGE("glimmer_sourcemod_dir: " ${CISM_SOURCEMOD_DIR})
# MESSAGE("Fortran Mod Files: ${SOURCEMODFILES}")

FOREACH( MODFILE ${SOURCEMODFILES})
  STRING(FIND ${MODFILE} / index REVERSE)
  MATH(EXPR index ${index}+1)
  STRING(SUBSTRING ${MODFILE} ${index} -1 filename)
  FOREACH( SOURCEFILE ${FORTRANSOURCES})
    STRING(REGEX MATCH ${filename} match_found ${SOURCEFILE})

    IF(match_found)
      MESSAGE("--SourceMod: removing ${SOURCEFILE} in favor of ${MODFILE}")
      LIST(REMOVE_ITEM FORTRANSOURCES ${SOURCEFILE})
   ENDIF()
  ENDFOREACH()
ENDFOREACH()
###

#IF (CISM_COUPLED)
  # enable removal of the autogenerated source files, when 'make clean' is done:
  # commented out for now, since it may create more problems than it solves when building
  # SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${FORTRAN_AUTOGEN_SOURCES}") 
#ENDIF()


SET(CMAKE_INCLUDE_PATH ${CISM_BINARY_DIR}/include)


#message("Fortransource ${FORTRANSOURCES}")

# include-dirs: Binary: all the .mod files; SOURCE: config.inc; 
#               libglide: glide_mask.inc

include_directories (${CISM_BINARY_DIR} ${CISM_SOURCE_DIR}
                     ${CISM_SOURCE_DIR}/libglide)


SET(LIBRARY_OUTPUT_PATH ${CISM_BINARY_DIR}/lib)

### Determine which fortran source files are fixed-form,
### which may require special compilation flags
FOREACH (SOURCEFILE ${FORTRANSOURCES} ${SOURCEMODFILES})
  STRING(REGEX MATCH "\\.[fF]$" match_found ${SOURCEFILE})
  IF(match_found)
    LIST(APPEND FIXEDSOURCES ${SOURCEFILE})
  ENDIF()
ENDFOREACH()
SET_SOURCE_FILES_PROPERTIES(${FIXEDSOURCES} PROPERTIES Fortran_FORMAT FIXED)

add_library(glimmercismfortran ${FORTRANSOURCES} ${FORTRAN_AUTOGEN_SOURCES} ${SOURCEMODFILES})

#####Build cism_driver executable (unless turned off) #####
OPTION(CISM_NO_EXECUTABLE "Set to  ON  to just build libraries (default:OFF)"  OFF)

IF (CISM_USE_CISM_FRONT_END)
   add_subdirectory(cism_driver)
ENDIF()


