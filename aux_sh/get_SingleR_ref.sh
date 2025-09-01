#! /usr/bin/env bash


# STAGE 0: SET LOCAL REFERENCE

#SBATCH -J get_SingleR_ref.sh
#SBATCH --cpus-per-task=48
#SBATCH --mem='1000gb'
#SBATCH --constraint=cal
#SBATCH --time=1-23:59:59
#SBATCH --error=job.ref.%J.err
#SBATCH --output=job.ref.%J.out
aux_opt=$1
if [ -z $SLURM_CPUS_PER_TASK ]; then
    CPU=1
else
    CPU=$SLURM_CPUS_PER_TASK
fi

source ~soft_bio_267/initializes/init_R

get_SingleR_ref.R --reference "$SingleR_ref" --version "$ref_version" \
                  --verbose $verbose --output "$ref_origin" \
                  --quiet $quiet --database "$ref_to_process" \
                  --ref_label "$ref_label" --cpu $CPU $aux_opt
