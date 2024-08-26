#!/bin/bash
for T in $(seq 10 30 120)
    do
    for theta in $(seq 270 30 360)
        do
            echo "$T $theta"
            ./a.out $T $theta $RANDOM
    done
done
