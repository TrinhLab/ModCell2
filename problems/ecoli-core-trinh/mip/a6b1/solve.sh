#!/bin/bash
PARENT=$(dirname `pwd`) 
mat2dat "${PARENT}/mip-lsgcp.mat" -o "mip-lsgcp.dat" -a 6 -b 1 -w ETH-100,BUTBUT-10 --free_w 0
runmc mip-lsgcp enumerate
# Cleanup unneccessary files
rm *.lp
rm *.dat
rm *.sol