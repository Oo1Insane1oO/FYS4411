# Project 1 - Hartree-Fock Methods

Uses Eigen class for linear algebra(finding eigenvalues), see https://www.eigen.tuxfamily.org

All files needed for running the Hartree-Fock program is given under src.

Data files(text files) for runs with different frequency values and different number of single particle orbitals and number of electrons are given under data. Sub-directories are for a given oscillator frequency while for instance the name w01_E5_e2.txt indicates a run with frequency 0.1, 5 single particle orbitals and 2 electrons.

## How to compile:
    1. make main

Report is given in text directory.

Python script for generating tables in latex-syntax is given as textable.py. A shell-script for running the python script for a given frequency is given as makeTables.sh.

A shell script for running for several of the mentioned parameters is given in run.sh.
