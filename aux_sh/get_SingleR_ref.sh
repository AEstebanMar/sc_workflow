#! /usr/bin/env bash


# STAGE 0: SET LOCAL REFERENCE

#SBATCH -J get_SingleR_ref.sh
#SBATCH --cpus-per-task=12
#SBATCH --mem='200gb'
#SBATCH --constraint=download
#SBATCH --time=1-23:59:59
#SBATCH --error=job.ref.%J.err
#SBATCH --output=job.ref.%J.out

aux_opt=$1
if [ -z $SLURM_CPUS_PER_TASK ]; then
    CPU=1
else
    CPU=$SLURM_CPUS_PER_TASK
fi

if [ "$aux_opt" == "--only_showcase" ]; then
    reference="$SingleR_ref"
    label="$ref_label"
else
    reference="$ref_to_process"
    label="$label_to_trim"
fi

source ~aestebanm/initializes/init_Hunter_dev
get_SingleR_ref.R --reference "$reference" --name "$ref_name" --cpu $CPU \
                  --version "$ref_version" --verbose $verbose --quiet $quiet \
                  --output "$ref_output_path" --ref_label "$label" $aux_opt
