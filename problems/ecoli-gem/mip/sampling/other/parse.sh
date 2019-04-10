# Keep only first reaction in each
tail -n +2 "fluxes.csv" | awk -F"," '{print $1}'  > rxnlist.txt
tail -n +2 "fluxes.csv" | awk -F"," '{print $1}' |sed 's/|/\n/g'  > rxnlistlong.txt