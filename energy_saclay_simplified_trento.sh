#!/bin/bash
export DYLD_LIBRARY_PATH=/usr/local/boost_1_67_0/stage/lib

cd /Users/Kianusch/Documents/Studium/Semester/WiSe1819/Bachelor-Arbeit/Heavy-Ion-Collision-Analysis

##### Generate Trento data

# Create the specific Trento script (name: PbPb.txt)
NEVENTS=100000
# DESTINATION=/Volumes/MAC/Trento/PbPb
DESTINATION=Trento/PbPb$NEVENTS
g++ -std=c++11 -O2 -o trento_generator trento_generator.cpp
./trento_generator $NEVENTS 0 20 $DESTINATION PbPb.txt

# Launch Trento
/Users/Kianusch/.local/bin/trento -c PbPb.txt

# Delete Trento script
rm -rf PbPb.txt




##### Compute magma data

cd /Users/Kianusch/Documents/Studium/Semester/WiSe1819/Bachelor-Arbeit/Heavy-Ion-Collision-Analysis/Saclay_simplified

# Compute Two-mode fluctuations from the simplified Saclay model
g++ -std=c++11 -O2 -lgsl -o SaclaySimpleTrento saclay_simplified_trento.cpp


centrality_min=20
centrality_max=21

m=0.14
Q0=1.24


./SaclaySimpleTrento $centrality_min $centrality_max $m $Q0 "../"$DESTINATION .dat $NEVENTS

# Plot result
#python3 plot_two_point_random_connected_fixed_m.py






