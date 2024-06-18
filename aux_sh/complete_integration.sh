#! /usr/bin/env bash


# STAGE 3 FULL INTEGRATION

#SBATCH -J complete_integration.sh
#SBATCH --cpus-per-task=12
#SBATCH --mem='200gb'
#SBATCH --constraint=cal
#SBATCH --time=0-23:00:00
#SBATCH --error=job.comp.%J.err
#SBATCH --output=job.comp.%J.out

