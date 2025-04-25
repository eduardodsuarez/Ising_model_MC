#!/bin/bash

## para compilador da GNU
g++ -c ising.cpp
g++ ising.o -o ising
## para compiladores intel
#icpc  -c  ising.cpp
#icpc ising.o -o ising_intel
