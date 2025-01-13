#!/bin/bash

# give execute permission to your script using chmod +x yourscript.sh before 
# running it with ./yourscript.sh


for i in {1..4} 
#for i in {1..2} 
#for i in {3..4} s
do
   julia save_expected_sca_lp.jl $i &
done

wait # Wait for all background processes to finish
echo "All processes 1-4 in save_expected_sca_lp.jl have completed."