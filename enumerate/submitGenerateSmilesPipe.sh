#!/bin/bash

source /lustre/fs6/lyu_lab/scratch/blam/soft/miniconda3/etc/profile.d/conda.sh
conda deactivate
conda activate analysis_env

function exists {
        env_name=$1
        desc=$2
        if [ -z "${!env_name}" ]; then
                echo "expected input arg: $env_name"
                echo "arg description: $desc" 
                failed=1
        fi
}

export WORKDIR=$PWD
exists TRIM_LIB
exists RXN_SCHEME
exists CHUNK
export ENUMERATE_BASE=/lustre/fs6/lyu_lab/scratch/blam/work/synthon/submit_repacking_mod

#Create smile combinations that are saved into json_chunks directory
python $ENUMERATE_BASE/combo_chunker.py $TRIM_LIB $CHUNK json_chunks

#Find the number of json files created for N array jobs
num_json=$(ls json_chunks/*json | wc -l)

mkdir log
LOGGING=$WORKDIR/log

#Output directory for process-json.sh
mkdir enumerated_smiles

sbatch $SBATCH_ARGS -o $LOGGING/%a.out -e $LOGGING/%a.er --array=1-"$num_json"%100 "$ENUMERATE_BASE"/process-json.sh -J enum_%a

echo "Finished generating smiles"


