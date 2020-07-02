#! /bin/bash

for i in $(seq 1 1 35); do 
    for alpha in 1 6; do
        echo $i $alpha
        python evolve.py `echo $i/100 | bc -l` $alpha
        python test.py mutated_a$alpha\_p$(printf "%03d" $i)
    done     
done
