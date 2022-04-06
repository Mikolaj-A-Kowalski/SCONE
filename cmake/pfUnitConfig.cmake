if(NOT BUILD_TESTS)
  message(FATAL_ERROR "This module should not be included if tests are not to be build.")
endif()

# FIND PYTHON INTERPRETER
# Sets PYTHONINTERP_FOUND & PYTHON_EXECUTABLE
find_package(PythonInterp REQUIRED)

#
# We need to set all variables provided by the find package
# See `FindPFUNIT.cmake`
#
if(EXTERNAL_PFUNIT)
  find_package(PFUNIT REQUIRED)
  add_library(pFUnit STATIC IMPORTED)
  set_property(TARGET pFUnit PROPERTY IMPORTED_LOCATION ${PFUNIT_LIBRARIES})

else()
  #
  # Using FetchContent allows us to make CMAKE download and build pFUnit by itself.
  # However, since we are using an old version of pFUnit we require some tricks:
  #  - We had to copy FindOpenMP_Fortran.cmake to cmake directory so it is in CMAKE_MODULE_PATH
  #  - We need to set key paths by hand
  #  - We need to tinker with the INCLUDE_DIRECTORIES of pfunit target so they do not contain
  #    missing directories
  #
  message(STATUS "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>")
  message(STATUS "Fetching pFUnit from: https://github.com/Goddard-Fortran-Ecosystem/pFUnit")
  message(STATUS "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>")

  # NOTE: Name MUST be in lowercase only. Otherwise FetchContent will not work withou any warnings
  FetchContent_Declare(
    pfunit_lib
    GIT_REPOSITORY  https://github.com/Goddard-Fortran-Ecosystem/pFUnit
    GIT_TAG         8d32a2cc63863588d55940bd5236ce68b0dc54c2 # v3.3.0 - said to be faster if specified by HASH (?) 
  )

  FetchContent_GetProperties(pfunit_lib)

  if (NOT pfunit_lib_POPULATED)
    FetchContent_Populate(pfunit_lib)

    # We need to EXCLUDE_FROM_ALL in order to build only the targets we need
    # This allows us to build pFUnit tests for example
    add_subdirectory(${pfunit_lib_SOURCE_DIR} ${pfunit_lib_BINARY_DIR} EXCLUDE_FROM_ALL)

    # We need to set approperiate locations by hand
    # May be a bit unstable, but we don't see a way around it at the moment
    set(PFUNIT_FOUND TRUE)
    set(PFUNIT_LIBRARIES ${pfunit_lib_BINARY_DIR}/source/libpfunit.a)
    set(PFUNIT_INCLUDE_DIRS ${pfunit_lib_SOURCE_DIR}/include)
    set(PFUNIT_MOD ${CMAKE_BINARY_DIR}/mod)
    set(PFUNIT_PREPROC ${pfunit_lib_SOURCE_DIR}/bin)

    # pfUnit tries to include '${CMAKE_BINARY_DIR}/source', which does not exist
    # we need to replace it with correct path '${pfunit_lib_BINARY_DIR}/source'
    get_target_property(pfunit_include pfunit INCLUDE_DIRECTORIES)
    list(REMOVE_ITEM pfunit_include "${CMAKE_BINARY_DIR}/source")
    list(APPEND pfunit_include "${pfunit_lib_BINARY_DIR}/source")
    set_property(TARGET pfunit PROPERTY INCLUDE_DIRECTORIES ${pfunit_include})

  endif()
endif()
