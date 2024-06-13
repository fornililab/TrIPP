#!/bin/bash 

#Run TrIPP on three replicas 
echo "Running pka calculations for three MD replicas" 
python test_installation.py 

#Check if there differences between files generated locally and model files 
pka_test=0
surf_test=0
for run in 1 2 3 
do 
echo "Testing pKa calculation for MD$run" 
output=$(diff 1AKI_MD$run'_comparison_file_pka.csv' 1AKI_MD$run'_local_file_pka.csv') 
if [ -n "$output" ]; then
  echo "WARNING: Detected differences between local calculations and comparison files." 
  echo "Please check the dependecies versions." 
else
  echo "Test successful"
  pka_test=$(( $pka_test + 1 ))
fi

echo "Testing surface exposure calculation for MD$run" 
output=$(diff 1AKI_MD$run'_comparison_file_surf.csv' 1AKI_MD$run'_local_file_surf.csv') 
if [ -n "$output" ]; then
  echo "WARNING: Detected differences between local calculations and comparison files." 
  echo "Please check the dependecies versions." 
else
  echo "Test successful"
  surf_test=$(( $surf_test + 1 ))
fi
done 

if [ $pka_test -eq 3 ]; then 
  echo "All pka calculations were successful" 
fi 

if [ $surf_test -eq 3 ]; then 
  echo "All surface exposure calculations were successful" 
fi 

if [ $pka_test -eq 3 -a $surf_test -eq 3 ]; then
  echo "All pka calculations and surface exposure calculations were successful"
  echo "Installation was successful" 
fi