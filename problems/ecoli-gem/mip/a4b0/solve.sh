#!/bin/bash
PARENT=$(dirname `pwd`)
PROBLEM="r-wgcp-17-tight"
TIMELIM=36000
mat2dat "${PARENT}/${PROBLEM}.mat" -o "${PROBLEM}.dat" --alpha 4 --free_w 0
runmc ${PROBLEM} benders $TIMELIM
# Cleanup unneccessary files
rm *.lp
rm *.dat
rm *.sol
