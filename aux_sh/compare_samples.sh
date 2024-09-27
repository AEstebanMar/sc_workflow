#! /usr/bin/env bash


# STAGE 2 SAMPLES COMPARISON

#SBATCH -J compare_samples.sh
#SBATCH --cpus-per-task=12
#SBATCH --mem='1000gb'
#SBATCH --constraint=cal
#SBATCH --time=2-23:59:59
#SBATCH --error=job.comp.%J.err
#SBATCH --output=job.comp.%J.out

# Setup

if [ -z $SLURM_CPUS_PER_TASK ]; then
    SLURM_CPUS_PER_TASK=1
fi

. ~soft_bio_267/initializes/init_ruby
. ~soft_bio_267/initializes/init_R
hostname

mkdir -p $output"/report"

if [ "$imported_counts" == "" ]; then
    cat $FULL_RESULTS/*/metrics > $experiment_folder'/metrics'
    cat $FULL_RESULTS/*/cellranger_metrics > $experiment_folder'/cellranger_metrics'
    create_metric_table.rb $experiment_folder'/metrics' sample $experiment_folder'/metric_table'
    create_metric_table.rb $experiment_folder'/cellranger_metrics' sample $experiment_folder'/cellranger_metric_table'

    compare_samples.R -o $output"/report" \
                      -m $experiment_folder'/metric_table' \
                      -l $experiment_folder'/metrics' \
                      -e $experiment_name \
                      --cellranger_metrics $experiment_folder'/cellranger_metric_table' \
                      --cellranger_long_metrics $experiment_folder'/cellranger_metrics'

fi

integration.R --output $output \
              --project_name $experiment_name \
              --filter $preproc_filter \
              --mincells $preproc_init_min_cells \
              --minfeats $preproc_init_min_feats \
              --minqcfeats $preproc_qc_min_feats \
              --percentmt $preproc_max_percent_mt \
              --normalmethod $preproc_norm_method \
              --scalefactor $preproc_scale_factor \
              --hvgs $preproc_select_hvgs \
              --ndims $preproc_pca_n_dims \
              --dimheatmapcells $preproc_pca_n_cells \
              --resolution $preproc_resolution \
              --exp_design $exp_design \
              --int_columns "$int_columns" \
              --count_path $FULL_RESULTS"/*/cellranger_0000/*" \
              --suffix "outs/filtered_feature_bc_matrix" \
              --cluster_annotation "$cluster_annotation" \
              --target_genes $target_genes \
              --cpu $SLURM_CPUS_PER_TASK \
              --imported_counts "$imported_counts" \
              --DEG_columns "$DEG_columns" \
              --cell_annotation "$cell_annotation" \
              --p_adj_cutoff $p_adj_cutoff \
              --verbose $verbose \
              --reduce $reduce \
              --samples_to_integrate $samples_to_process
