#!/bin/bash

# give execute permission to your script using chmod +x yourscript.sh before 
# running it with ./yourscript.sh


#for lowpass = true only
# run 1-10 with 100 each
for i in {1..10} 
#for i in {1..2} 
#for i in {3..4} 
#for i in {5..6} 
#for i in {7..10} 
do
   julia save_expected_hurst_lp.jl $i &
done

wait # Wait for all background processes to finish
echo "All processes 1-10 in save_expected_hurst_lp.jl have completed."