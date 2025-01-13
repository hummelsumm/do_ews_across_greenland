#!/bin/bash

# give execute permission to your script using chmod +x yourscript.sh before 
# running it with ./yourscript.sh

for i in {1..4}
do
   julia save_20y_wavelet.jl $i &
done

wait # Wait for all background processes to finish
echo "All 4 processes in save_20y_wavelet.jl have completed."