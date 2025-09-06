
sbatch <<EOF
#!/bin/bash
#SBATCH --partition=lyu_docking,lyu
#SBATCH --output=./output.log
#SBATCH --error=./myjob.err
#SBATCH --job-name=genSmiles


export TRIM_LIB=$1
export RXN_SCHEME=$2
export CHUNK=$3
export SBATCH_ARGS="--partition=lyu_docking,lyu"

sh /lustre/fs6/lyu_lab/scratch/blam/work/synthon/submit_repacking_mod/submitGenerateSmilesPipe.sh
EOF


