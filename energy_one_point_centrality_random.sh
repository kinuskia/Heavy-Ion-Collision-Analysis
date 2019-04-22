#!/bin/bash
export DYLD_LIBRARY_PATH=/usr/local/boost_1_67_0/stage/lib
cd /Users/Kianusch/Documents/Studium/Semester/WiSe1819/Bachelor-Arbeit/Heavy-Ion-Collision-Analysis


# Create the specific Trento script (name: PbPb.txt)
NEVENTS=10000
# DESTINATION=/Volumes/MAC/Trento/PbPb
DESTINATION=Trento/PbPb$NEVENTS
g++ -std=c++11 -O2 -o trento_generator trento_generator.cpp
./trento_generator $NEVENTS 0 20 $DESTINATION PbPb.txt



# Launch Trento
/Users/Kianusch/.local/bin/trento -c PbPb.txt

# Delete Trento script
rm -rf PbPb.txt
	
# Evaluate centrality of Trento files with name Trento/PbPb.dat, 
g++ -std=c++11 -lgsl -O2 -o OnePointCentralityRandom analysis_one_point_centrality_random.cpp
./OnePointCentralityRandom $DESTINATION .dat $NEVENTS


