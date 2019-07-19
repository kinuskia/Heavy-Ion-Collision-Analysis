#!/bin/bash
export DYLD_LIBRARY_PATH=/usr/local/boost_1_67_0/stage/lib
cd /Users/Kianusch/Documents/Studium/Semester/WiSe1819/Bachelor-Arbeit/Heavy-Ion-Collision-Analysis

NEVENTS=1000
DESTINATION=Trento/PbPb

for b in $(seq -w 0 0) # for each impact parameter
do
	# Create the specific Trento script (name: PbPb.txt)
	./trento_generator $NEVENTS $b $b $DESTINATION$b PbPb.txt

	# Launch Trento
	/Users/Kianusch/.local/bin/trento -c PbPb.txt

	# Delete Trento script
	rm -rf PbPb.txt
	
	# Do two-point functions on Trento files with name Trento+impact_parameter+.dat, 
	# each folder having 1000 files
	for m in $(seq -w 0 9) # for each parameter m
	do
	./Two_Point_Correlators $DESTINATION$b .dat $b $m $NEVENTS
	done
done