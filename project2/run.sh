#!/bin/bash
for E in 2 6 12 20
do
    for OMEGA in 1 5 10
    do
        mpirun -np 8 ./vmc_main `bc<<<"scale=2;"$OMEGA"/10"` $E 1000000 0.02 0 1 1 1 "data/w0"$OMEGA"/N"$E"/w"$OMEGA"_N"$E"_"
    done
done

for E in 2 6 12 20
do
    for OMEGA in 1 5
    do
        mpirun -np 8 ./vmc_main `bc<<<"scale=2;"$OMEGA"/100"` $E 1000000 0.02 0 1 1 1 "data/w00"$OMEGA"/N"$E"/w"$OMEGA"_N"$E"_"
    done
done
