#!/bin/bash
for OMEGA in 1 2 3 4 5
do
    for E in 2 6 12 20
    do
        mpirun -np 4 ./vmc_main `bc<<<"scale=2;"$OMEGA"/10"` $E 1000000 0.00035 0 1 1 1 "data/w0"$OMEGA"/N"$E"/w"$OMEGA"_N"$E"_"
    done
done
