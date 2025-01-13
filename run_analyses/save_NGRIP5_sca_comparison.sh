#!/bin/bash

# give execute permission to your script using chmod +x yourscript.sh before 
# running it with ./yourscript.sh

for i in {1..9} 
#for i in {1..3}
#for i in {4..6}
#for i in {7..9}
do
   julia save_NGRIP5_sca_comparison.jl $i &
done

wait # Wait for all background processes to finish
echo "All 9 processes in save_NGRIP5_sca_comparison.jl have completed."