#!/bin/bash
export DYLD_LIBRARY_PATH=/usr/local/boost_1_67_0/stage/lib
cd /Users/Kianusch/Documents/Studium/Semester/WiSe1819/Bachelor-Arbeit/Heavy-Ion-Collision-Analysis/Saclay_simplified


# Compute Two-mode fluctuations from the simplified Saclay model
g++ -std=c++11 -O2 -lgsl -o SaclaySimple saclay_simplified.cpp

# m=0.6
# Q0=1.1

# m=0.15
# Q0=2.4

centrality_min=20
centrality_max=21
m=0.15
Q0=1.1


./SaclaySimple $centrality_min $centrality_max $m $Q0

# Plot result
#python3 plot_two_point_random_connected_fixed_m.py






