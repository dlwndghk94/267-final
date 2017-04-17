#!/bin/bash

echo "Running standard algorithm \n"
rm output.csv
rm *.out
g++ particle_filter.cpp
./a.out

echo "Running omp algorithm"
rm output.csv
rm *.out
export OMP_NUM_THREADS=6
g++ -fopenmp particle_filter_omp.cpp
./a.out

# cd visualizer
# python particle_visualizer.py '../output.csv' 100 100
