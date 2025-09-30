#! /usr/bin/env bash

#SBATCH -J sketch_sc.sh
#SBATCH --cpus-per-task=24
#SBATCH --mem='400gb'
#SBATCH --constraint=cal
#SBATCH --time=1-00:00:00
#SBATCH --error=job.sing.%J.err
#SBATCH --output=job.sing.%J.out


if [ "$1" == "" ]; then
	echo "ERROR: UNDEFINED SCRIPT"
	exit 1
fi

$CODE_PATH/singularity/seurat_image.sh exec $1 $SLURM_CPUS_PER_TASK
