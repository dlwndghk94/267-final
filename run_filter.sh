#!/bin/bash


echo "Running standard algorithm..."
rm output.csv
rm particle_filter_timed.out
g++ -O3 -o particle_filter_timed.out particle_filter_timed.cpp
#./a.out

echo "Running naive algorithm..."
rm output.csv
rm particle_filter_serial_naive.out
g++ -O3 -o particle_filter_serial_naive.out particle_filter_serial_naive.cpp
./a.out

#echo "\nRunning omp algorithm, 4 threads"
#rm output.csv
#rm *.out
#export OMP_NUM_THREADS=4
#g++ -fopenmp -O3 particle_filter_omp.cpp
#./a.out

#echo "\nRunning omp algorithm, 8 threads"
#rm output.csv
#rm *.out
#export OMP_NUM_THREADS=8
#g++ -fopenmp -O3 particle_filter_omp.cpp
#./a.out


# echo "\nRunning omp algorithm, 4 threads"
# rm output.csv
# rm *.out
# export OMP_NUM_THREADS=4
# g++ -fopenmp -O3 particle_filter_omp.cpp
# ./a.out

# echo "\nRunning omp algorithm, 8 threads"
# rm output.csv
# rm *.out
# export OMP_NUM_THREADS=8
# g++ -fopenmp -O3 particle_filter_omp.cpp
# ./a.out

# echo "\nRunning omp algorithm, 8 threads"
# rm output.csv
# rm *.out
# export OMP_NUM_THREADS=8
# g++ -fopenmp particle_filter_omp.cpp
# ./a.out

# cd visualizer
# python particle_visualizer.py '../output.csv' 100 100
