#!/bin/bash

# give execute permission to your script using chmod +x yourscript.sh before 
# running it with ./yourscript.sh


for i in {1..4} 
do
   julia save_expected_csd.jl $i &
done

wait # Wait for all background processes to finish
echo "All processes have completed."