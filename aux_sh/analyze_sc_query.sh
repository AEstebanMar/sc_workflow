#! /usr/bin/env bash

#! /usr/bin/env bash


# STAGE 5 QUERY GENE ANALYSIS

#SBATCH -J analyze_sc_query.sh
#SBATCH --cpus-per-task=12
#SBATCH --mem='300gb'
#SBATCH --constraint=cal
#SBATCH --time=0-5:00:00
#SBATCH --error=job.query.%J.err
#SBATCH --output=job.query.%J.out

if [ -z $SLURM_CPUS_PER_TASK ]; then
    SLURM_CPUS_PER_TASK=1
fi

. ~aestebanm/software/inits/init_Hunter_dev
hostname
mkdir -p $output"/report"
analyze_sc_query.R --input $output \
                   --extra_columns "$extra_columns" \
            	   --target_genes $target_genes \
            	   --output $output \
                   --sigfig $sigfig \
            	   --name $experiment_name \
            	   --verbose $verbose \
            	   --cpu $SLURM_CPUS_PER_TASK \
                   --min_counts $min_counts


