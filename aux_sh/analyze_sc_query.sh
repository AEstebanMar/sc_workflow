#! /usr/bin/env bash

#! /usr/bin/env bash


# STAGE 5 QUERY GENE ANALYSIS

#SBATCH -J analyze_sc_query.sh
#SBATCH --cpus-per-task=12
#SBATCH --mem='300gb'
#SBATCH --constraint=cal
#SBATCH --time=0-5:00:00
#SBATCH --error=job.DEG.%J.err
#SBATCH --output=job.DEG.%J.out

if [ -z $SLURM_CPUS_PER_TASK ]; then
    SLURM_CPUS_PER_TASK=1
fi

. ~aestebanm/initializes/init_Hunter_dev
hostname
mkdir -p $output"/report"

analyze_sc_query.R --DEG_target $DEG_target \
            	   --input $output/counts \
            	   --target_genes $target_genes \
            	   --p_val_cutoff $DEG_p_val_cutoff \
            	   --min_avg_log2FC $min_avg_log2FC \
            	   --min_cell_proportion $min_cell_proportion \
            	   --min_counts $min_counts \
            	   --output $output \
            	   --name $experiment_name \
            	   --verbose $verbose \
            	   --cpu $SLURM_CPUS_PER_TASK


