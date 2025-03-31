#! /usr/bin/env bash

#SBATCH -J seurat.sh
#SBATCH --cpus-per-task=12
#SBATCH --mem='700gb'
#SBATCH --constraint=cal
#SBATCH --time=0-20:00:00
#SBATCH --error=job.seu.%J.err
#SBATCH --output=job.seu.%J.out


if [ "$1" == "" ]; then
	echo "ERROR: UNDEFINED SCRIPT"
	exit 1
fi

./seurat.sh exec "$config_file 2"
