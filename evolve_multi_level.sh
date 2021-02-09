#! /bin/bash

alpha=22

for i in $(seq 1 2 70); do 
    #for alpha in 1 6 12 1000; do
        echo $i $alpha
        python evolve.py `echo $i/100 | bc -l` $alpha
        python compute_GND_AAD.py mutated_a$alpha\_aad$(printf "%03d" $i) 
    #done     
done
