#! /bin/bash

for i in $(seq 1 1 30); do 
    echo $i
    python evolve.py `echo $i/100 | bc -l`
    python test.py mutated_$(printf "%03d" $i)
done
