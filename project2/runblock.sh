#!/bin/bash
for E in 2 6 12 20
do
    for OMEGA in 1 5 10
    do
        python "scripts/readBinary.py" "data/w0"$OMEGA"/N"$E"/" 1000
    done
done

for E in 2 6 12 20
do
    for OMEGA in 1 5
    do
        python "scripts/readBinary.py" "data/w00"$OMEGA"/N"$E"/" 1000
    done
done
