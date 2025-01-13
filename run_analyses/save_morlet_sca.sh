#!/bin/bash

# give execute permission to your script using chmod +x yourscript.sh before 
# running it with ./yourscript.sh

for i in {1..6} 
#for i in {1..4}
#for i in {5..6}
do
   julia save_morlet_sca.jl $i &
done

wait # Wait for all background processes to finish
echo "All 6 processes in save_morlet_sca.jl have completed."