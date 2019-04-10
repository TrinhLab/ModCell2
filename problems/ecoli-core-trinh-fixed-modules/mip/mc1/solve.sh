#!/bin/bash
PARENT=$(dirname `pwd`) 
mat2dat "${PARENT}/mip-lsgcp.mat" -o "mip-lsgcp.dat" --alpha 5
runmc mip-lsgcp enumerate
# Cleanup unneccessary files
rm *.lp
rm *.dat
rm *.sol