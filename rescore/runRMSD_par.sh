#!/bin/bash

while IFS= read -r mol2; do
    # Submit the job to SLURM with the appropriate inputs
    sbatch --partition=hpc --wrap="python calculateRMSD.py \"$mol2\" trimmedSynthons.mol2"
done < flist
