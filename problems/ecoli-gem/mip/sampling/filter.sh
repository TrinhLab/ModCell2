#!/bin/bash
rfile="other/rxnlist.txt"
for sfile in /samples/*.csv; do
    sed 's/^/^/' "$rfile" | xargs -ix grep x "$sfile" >> /fsamples/$(basename "$filename" .csv).csv
done