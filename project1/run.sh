#!/bin/bash
for OMEGA in 1 2 3 4 5
do
    for SHELL in 5 6 7 8 9 10 11
    do
        for E in 2 6 12 20
        do
            ./main `bc<<<"scale=2;"$OMEGA"/10"` $SHELL $E 1000 0 >> "data/w0"$OMEGA"/w0"$OMEGA"_E"$SHELL"_e"$E".txt"
        done
    done
done
