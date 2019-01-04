#!/bin/bash
PARENT=$(dirname `pwd`) 
mat2dat "${PARENT}/mip-lsgcp.mat" -o "mip-lsgcp.dat" -a 5 -b 1
runmc mip-lsgcp enumerate
# Cleanup unneccessary files
rm *.lp
rm *.dat
rm *.sol