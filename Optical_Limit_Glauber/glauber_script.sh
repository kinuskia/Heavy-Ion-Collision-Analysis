#!/bin/bash
export DYLD_LIBRARY_PATH=/usr/local/boost_1_67_0/stage/lib
cd /Users/Kianusch/Documents/Studium/Semester/WiSe1819/Bachelor-Arbeit/Heavy-Ion-Collision-Analysis/Optical_Limit_Glauber/

FILENAME=mult-b_x

for x in $(seq -w 0 1 020) # value of x
do
	NEWFILENAME=$FILENAME$x
	./Glauber $x $NEWFILENAME
	
done
