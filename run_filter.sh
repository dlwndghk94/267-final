#!/bin/bash


#echo "Running standard algorithm..."
#rm output.csv
#rm *.out
#g++ particle_filter_timed.cpp
#./a.out

echo "\nRunning omp algorithm..."
rm output.csv
rm *.out
export OMP_NUM_THREADS=8
g++ -fopenmp particle_filter_omp.cpp
./a.out

#echo "\nRunning GPU algorithm..."
#rm output.csv
#rm *.out
#nvcc particle_filter_gpu.cu
#./a.out

# cd visualizer
# python particle_visualizer.py '../output.csv' 100 100
