#!/bin/bash
PARENT=$(dirname `pwd`)
echo 
mat2dat "${PARENT}/mip-lsgcp.mat" -o "mip-lsgcp.dat" --alpha 5 -w PRO-10,PROBUT-100
runmc mip-lsgcp enumerate
# Cleanup unneccessary files
rm *.lp
rm *.dat
rm *.sol