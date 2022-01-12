cmake_minimum_required (VERSION 3.5)
project (caobs_f2c_tests LANGUAGES C Fortran)
set(CMAKE_BUILD_TYPE Debug)

# Test 1 - Fortran version of refvap from caobs_main, and a simple
# test code.
add_executable(f_test_refvap src/fortran_tests/test_refvap.f90 src/caobs_main/refvap.f90)
# The following line should be changed and uncommented when a suitable
# C conversion is ready to be built and tested.
#add_executable(c_test_refvap src/c_tests/test_refvap.c src/caobs_main/refvap.c)
