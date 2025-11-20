#! /usr/bin/env bash


# STAGE 4 DEG ANALYSIS

#SBATCH -J sc_Hunter.sh
#SBATCH --cpus-per-task=36
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

sc_Hunter.R --targets_folder $TARGETS_FOLDER \
            --input $output/counts \
            --target_genes $target_genes \
            --p_val_cutoff $DEG_p_val_cutoff \
            --min_avg_log2FC $min_avg_log2fc \
            --min_cell_proportion $min_cell_proportion \
            --simple_DEG_pct $simple_DEG_pct \
            --min_counts $min_counts \
            --output $output \
            --name $experiment_name \
            --verbose $verbose \
            --cpu $SLURM_CPUS_PER_TASK \
            --top_N $top_N \
            --DE_method "$ref_de_method" \
            --minpack_common $minpack_common
