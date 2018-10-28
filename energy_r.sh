#!/bin/bash

for b in $(seq -w 0 19)
do
	./radial_analysis ../Trento-Example/PbPb .dat $b 1000
done
