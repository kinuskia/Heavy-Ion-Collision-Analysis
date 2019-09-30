#!/bin/bash
export DYLD_LIBRARY_PATH=/usr/local/boost_1_67_0/stage/lib
cd /Users/Kianusch/Documents/Studium/Semester/WiSe1819/Bachelor-Arbeit/Heavy-Ion-Collision-Analysis/Saclay


# Compute Two-mode fluctuations from the simplified Saclay model
g++ -std=c++11 -O2 -lgsl -o Saclay saclay.cpp


#centrality_min=20
#centrality_max=21
m=1.4e-1 #1e-7 .. 1e-2, 2e-2 2.05e-2 2.09e-2
#n_grid=10

for centrality_min in $(seq 0 99)
do
	centrality_max=$((centrality_min+1))
	centrality=$centrality_min"-"$centrality_max
	for n_grid in 41
	do
		for n_azim in 64
		do
				# create necessary directories
				destination=output/$centrality/Nr$n_grid/Nm$n_azim/m$m
				mkdir -p $destination
				./Saclay $centrality_min $centrality_max $n_grid $n_azim $destination $m
		done
	done
done


# Plot result
#python3 plot_two_point_random_connected_fixed_m.py






