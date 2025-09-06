#!/bin/bash

WORKDIR=$PWD

for id in `cat rxn.ids`; do
    cd $WORKDIR

    rm -r "$id"

    echo "$id"
    mkdir $id
    grep $id 2023-04_REAL_synthons_FIXED.txt > "$id"/"$id"_synthons.txt
    grep $id 2023-04_REAL_reactions_SMARTS.txt > "$id"/"$id"_rxn.txt
    cp generateSmiles.sh "$id"

    cd $id
    
    bash generateSmiles.sh "$id"_synthons.txt "$id"_rxn.txt 2000000

done
