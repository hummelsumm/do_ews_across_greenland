#!/bin/bash

# give execute permission to your script using chmod +x yourscript.sh before 
# running it with ./yourscript.sh

for i in {1..7} 
#for i in {1..4}
#for i in {5..7}
do
   julia save_NGRIP5_csd.jl $i &
done

wait # Wait for all background processes to finish
echo "All 7 processes in save_NGRIP5_csd.jl have completed."