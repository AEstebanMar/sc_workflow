#! /usr/bin/env bash

#SBATCH -J seurat.sh
#SBATCH --cpus-per-task=64
#SBATCH --mem='1000gb'
#SBATCH --constraint=cal
#SBATCH --time=4-22:00:00
#SBATCH --error=job.seu.%J.err
#SBATCH --output=job.seu.%J.out


if [ "$1" == "" ]; then
	echo "ERROR: UNDEFINED CONFIG FILE"
	exit 1
fi
export config_file=$1

./seurat.sh exec "../daemon.sh $config_file 2"
