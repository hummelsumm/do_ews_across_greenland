#!/bin/bash

# give execute permission to your script using chmod +x yourscript.sh before 
# running it with ./yourscript.sh

for i in {1..11} 
#for i in {1..4}
#for i in {5..8}
#for i in {9..11}
do
   julia save_NGRIP5_all_scales_hurst.jl $i &
done

wait # Wait for all background processes to finish
echo "All 11 processes in save_NGRIP5_all_scales_hurst.jl have completed."