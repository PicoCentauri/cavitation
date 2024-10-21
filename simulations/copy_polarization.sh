#!/bin/bash
echo $1

for RATE in $(find . -maxdepth 1 -mindepth 1 -executable); do
    mkdir -p $1/$RATE
    cp $RATE/grompp.mdp $1/$RATE
done

cp topol.top posreCH3.itp $1
