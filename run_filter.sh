#!/bin/bash

rm output.csv
g++ particle_filter.cpp
./a.out


# cd visualizer
# python particle_visualizer.py '../output.csv' 100 100

