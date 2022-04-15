#
# SETS TWO VARIABLES
#  SCONE_COMPILER_FLAGS -> FLAGS used when compiling SCONE library and executable
#  TESTS_COMPILER_FLAGS -> FLAGS used when compiling SCONE library and executable
#

# Make sure to clean the working variables
unset(SCONE_FLAGS_LOC CACHE)
unset(TESTS_FLAGS_LOC CACHE)

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(SCONE_FLAGS_LOC -std=f2008 -O3 -g -pedantic -Wall -Wno-unused-dummy-argument -cpp)
  # Tests do not require optimisation, we also need to disable warnings
  set(TESTS_FLAGS_LOC -std=f2008 -O0 -g -cpp)

  if(DEBUG)
    list(APPEND SCONE_FLAGS_LOC  -fcheck=bounds,do,mem,pointer -Waliasing)
    list(APPEND TESTS_FLAGS_LOC  -fcheck=bounds,do,mem,pointer -Waliasing)
  endif()

else ()
  message(STATUS ${CMAKE_Fortran_COMPILER_ID} EQUAL "GNU")
  message(FATAL_ERROR "Trying to use unsupported Fortran compiler:
  ${CMAKE_Fortran_COMPILER_ID}
  ${CMAKE_Fortran_COMPILER} ${CMAKE_Fortran_COMPILER_VERSION}" )

endif()

# EXPORT GLOBAL Variables
set(SCONE_COMPILER_FLAGS ${SCONE_FLAGS_LOC})
set(TESTS_COMPILER_FLAGS ${TESTS_FLAGS_LOC})
