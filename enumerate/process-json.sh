#!/bin/bash

ID=$(($SLURM_ARRAY_TASK_ID - 1))
json_file="$PWD"/json_chunks/*-"$ID".json

echo Generating smiles file for $json_file
python $ENUMERATE_BASE/build_combos.py $json_file $TRIM_LIB $RXN_SCHEME enumerated_smiles
echo Finished!
