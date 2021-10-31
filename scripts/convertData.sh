#!/bin/bash

for line in $(cat ./datasets/datasets.txt)
do
    echo $"Input the node label number for dataset $line:"
    read x
    echo $"Input the edge label number for dataset $line:"
    read y
    echo $"process dataset $line..."
    python convertDataToOurFormatSx.py $line $x $y
    echo -e $"dataset $line processed.\n"
    sleep 1
done
