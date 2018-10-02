#!/bin/sh

# build fftMPI library

echo "---------------------------"
echo "Building fftmpi library ..."
echo "---------------------------"
echo

cd ../src
make
cd ../test

# build test apps

echo
echo "---------------------------"
echo "Building 6 test apps ..."
echo "---------------------------"
echo

make

# run test apps in each language

echo
echo "---------------------------"
echo "Running C++ test apps in serial ..."
echo "---------------------------"
echo

test3d -g 64 64 64 -n 50
test2d -g 128 128 -n 100

echo
echo "---------------------------"
echo "Running C test apps in serial ..."
echo "---------------------------"
echo

test3d_c -g 64 64 64 -n 50
test2d_c -g 128 128 -n 100

echo
echo "-----------------------"
echo "Running Fortran test apps in serial ..."
echo "-----------------------"
echo

test3d_f90 -g 64 64 64 -n 50
test2d_f90 -g 128 128 -n 100

echo
echo "-----------------------"
echo "Running Python test apps in serial ..."
echo "-----------------------"
echo

python test3d.py -g 64 64 64 -n 50
python test2d.py -g 128 128 -n 100

echo
echo "---------------------------"
echo "Running C++ test apps in parallel ..."
echo "---------------------------"
echo

mpirun -np 8 test3d -g 64 64 64 -n 50
mpirun -np 8 test2d -g 128 128 -n 100

echo
echo "---------------------------"
echo "Running C test apps in parallel ..."
echo "---------------------------"
echo

mpirun -np 8 test3d_c -g 64 64 64 -n 50
mpirun -np 8 test2d_c -g 128 128 -n 100

echo
echo "-----------------------"
echo "Running Fortran test apps in parallel ..."
echo "-----------------------"
echo

mpirun -np 8 test3d_f90 -g 64 64 64 -n 50
mpirun -np 8 test2d_f90 -g 128 128 -n 100

echo
echo "-----------------------"
echo "Running Python test apps in parallel ..."
echo "-----------------------"
echo

mpirun -np 8 python test3d.py -g 64 64 64 -n 50
mpirun -np 8 python test2d.py -g 128 128 -n 100
