#!/bin/bash
export DYLD_LIBRARY_PATH=/usr/local/boost_1_67_0/stage/lib
cd /Users/Kianusch/Documents/Studium/Semester/WiSe1819/Bachelor-Arbeit/Heavy-Ion-Collision-Analysis/Saclay_simplified


# Compute Two-mode fluctuations from the simplified Saclay model
g++ -std=c++11 -O2 -lgsl -o SaclaySimple saclay_simplified.cpp
./SaclaySimple 

# Plot result
python3 plot_two_point_random_connected_fixed_m.py


