cmake_minimum_required (VERSION 3.5)
project (caobs_f2c_tests LANGUAGES C Fortran)
set(CMAKE_BUILD_TYPE Debug)

add_executable(f_test_refvap src/fortran_tests/test_refvap.f90 src/caobs_main/refvap.f90)
