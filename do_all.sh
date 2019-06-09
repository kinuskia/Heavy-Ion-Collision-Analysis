#!/bin/bash
export DYLD_LIBRARY_PATH=/usr/local/boost_1_67_0/stage/lib
cd /Users/Kianusch/Documents/Studium/Semester/WiSe1819/Bachelor-Arbeit/Heavy-Ion-Collision-Analysis

./energy_one_point_centrality.sh

./energy_two_point_centrality_random_fixed_m.sh

python3 plot_one_point_centrality.py

python3 plot_two_point_centrality_random_fixed_m.py


