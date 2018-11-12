#!/bin/bash
export DYLD_LIBRARY_PATH=/usr/local/boost_1_67_0/stage/lib
cd /Users/Kianusch/Documents/Studium/Semester/WiSe1819/Bachelor-Arbeit/Heavy-Ion-Collision-Analysis

for b in $(seq -w 0 11) # for each impact parameter
do
	# Create the specific Trento script (name: PbPb.txt)
	./trento_generator $b PbPb.txt

	# Launch Trento
	/Users/Kianusch/.local/bin/trento -c PbPb.txt

	# Delete Trento script
	rm -rf PbPb.txt
	
	# Do radial analysis on Trento files with name Trento+impact_parameter+.dat, 
	# each folder having 1000 files
	./radial_analysis Trento/PbPb .dat $b 1000
done
