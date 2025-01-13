#!/bin/bash

# give execute permission to your script using chmod +x yourscript.sh before 
# running it with ./yourscript.sh

# for lowpass = false
for i in {1..13} 
#for i in {1..4} 
#for i in {5..8} 
#for i in {9..12} 
#for i in {13..13} 
do
   julia save_expected_hurst_no_lp.jl $i &
done

wait # Wait for all background processes to finish
echo "All processes 1-10 in save_expected_hurst_no_lp.jl have completed."
