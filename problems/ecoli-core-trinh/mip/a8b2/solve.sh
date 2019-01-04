#!/bin/bash
PARENT=$(dirname `pwd`) 
mat2dat "${PARENT}/mip-lsgcp.mat" -o "mip-lsgcp.dat" -a 8 -b 2 --free_w 0
runmc mip-lsgcp enumerate
# Cleanup unneccessary files
rm *.lp
rm *.dat
rm *.sol