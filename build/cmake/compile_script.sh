#!/bin/bash

# Script to compile the SUMMA model using cmake
# This will call the CMakeLists.txt file in one directory above the current
# directory, and create a directory cmake_build where all of the build files
# will be stored. This script will build Summa in parallel.

rm -rf cmake_build



#cmake -B cmake_build -S ../. -DUSE_MIZUROUTE=ON
cmake -B cmake_build -S ../. \
  -DUSE_MIZUROUTE=ON \
  -DNetCDF_C_INCLUDE_DIR=/usr/include \
  -DNetCDF_F90_INCLUDE_DIR=/usr/lib64/gfortran/modules
  
cmake --build cmake_build --target all -j 