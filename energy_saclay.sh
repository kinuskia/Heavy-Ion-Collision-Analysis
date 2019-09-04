#!/bin/bash
export DYLD_LIBRARY_PATH=/usr/local/boost_1_67_0/stage/lib
cd /Users/Kianusch/Documents/Studium/Semester/WiSe1819/Bachelor-Arbeit/Heavy-Ion-Collision-Analysis/Saclay


# Compute Two-mode fluctuations from the simplified Saclay model
g++ -std=c++11 -O2 -lgsl -o Saclay saclay.cpp


centrality_min=20
centrality_max=21
centrality=$centrality_min"-"$centrality_max
#n_grid=10

for n_grid in 29
do
	for n_azim in 64
	do
		# create necessary directories
		destination=output/$centrality/Nr$n_grid/Nm$n_azim
		mkdir -p $destination
		./Saclay $centrality_min $centrality_max $n_grid $n_azim $destination
	done
done



# Plot result
#python3 plot_two_point_random_connected_fixed_m.py






