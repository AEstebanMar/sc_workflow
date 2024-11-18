#! /usr/bin/env bash


# STAGE 0: SET LOCAL REFERENCE

#SBATCH -J get_SingleR_ref.sh
#SBATCH --cpus-per-task=24
#SBATCH --mem='1000gb'
#SBATCH --constraint=cal
#SBATCH --time=0-23:59:59
#SBATCH --error=job.comp.%J.err
#SBATCH --output=job.comp.%J.out
aux_opt=$1
if [ -z $SLURM_CPUS_PER_TASK ]; then
    CPU=1
else
    CPU=$SLURM_CPUS_PER_TASK
fi
get_SingleR_ref.R --reference "$SingleR_ref" --version "$ref_version" \
                  --output $refs_path --verbose $verbose \
                  --quiet $quiet --database "$ref_origin" \
                  --ref_label $ref_label --cpu $CPU $aux_opt
