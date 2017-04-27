#!/bin/bash


echo "Running standard algorithm..."
rm output.csv
rm *.out
g++ -O3 particle_filter_timed.cpp
./a.out

echo "\nRunning omp algorithm, 2 threads"
rm output.csv
rm *.out
export OMP_NUM_THREADS=2
g++ -fopenmp -O3 particle_filter_omp.cpp
./a.out

echo "\nRunning omp algorithm, 4 threads"
rm output.csv
rm *.out
export OMP_NUM_THREADS=4
g++ -fopenmp -O3 particle_filter_omp.cpp
./a.out

echo "\nRunning omp algorithm, 8 threads"
rm output.csv
rm *.out
export OMP_NUM_THREADS=8
g++ -fopenmp -O3 particle_filter_omp.cpp
./a.out

# echo "\nRunning omp algorithm, 8 threads"
# rm output.csv
# rm *.out
# export OMP_NUM_THREADS=8
# g++ -fopenmp particle_filter_omp.cpp
# ./a.out

# cd visualizer
# python particle_visualizer.py '../output.csv' 100 100
