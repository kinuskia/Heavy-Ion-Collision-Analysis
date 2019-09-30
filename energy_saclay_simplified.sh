#!/bin/bash
export DYLD_LIBRARY_PATH=/usr/local/boost_1_67_0/stage/lib
cd /Users/Kianusch/Documents/Studium/Semester/WiSe1819/Bachelor-Arbeit/Heavy-Ion-Collision-Analysis/Saclay_simplified


# Compute Two-mode fluctuations from the simplified Saclay model
g++ -std=c++11 -O2 -lgsl -o SaclaySimple saclay_simplified.cpp

# m=0.6
# Q0=1.1

# m=0.15
# Q0=2.4

#centrality_min=0
#centrality_max=1
# m=0.15
# Q0=1.1
m=0.14
Q0=1.24

for centrality_min in $(seq 0 99)
do
	centrality_max=$((centrality_min+1))
	centrality=$centrality_min"-"$centrality_max
	destination=output/$centrality
	mkdir -p $destination
	./SaclaySimple $centrality_min $centrality_max $m $Q0 $destination
done




# Plot result
#python3 plot_two_point_random_connected_fixed_m.py






