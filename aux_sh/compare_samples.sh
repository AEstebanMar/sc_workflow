#! /usr/bin/env bash

# Sergio Alias, 20230530
# Last modified 20230721

# STAGE 2 SAMPLES COMPARISON

#SBATCH -J compare_samples.sh
#SBATCH --cpus-per-task=12
#SBATCH --mem='200gb'
#SBATCH --constraint=cal
#SBATCH --time=0-23:00:00
#SBATCH --error=job.comp.%J.err
#SBATCH --output=job.comp.%J.out

# Setup

. ~soft_bio_267/initializes/init_ruby
. ~soft_bio_267/initializes/init_R
hostname

mkdir -p $report_folder


cat $FULL_RESULTS/*/metrics > $experiment_folder'/metrics'
cat $FULL_RESULTS/*/cellranger_metrics > $experiment_folder'/cellranger_metrics'
create_metric_table.rb $experiment_folder'/metrics' sample $experiment_folder'/metric_table'
create_metric_table.rb $experiment_folder'/cellranger_metrics' sample $experiment_folder'/cellranger_metric_table'

/usr/bin/time $CODE_PATH/scripts/compare_samples.R -o $report_folder \
                                                   -m $experiment_folder'/metric_table' \
                                                   -l $experiment_folder'/metrics' \
                                                   -e $experiment_name \
                                                   --cellranger_metrics $experiment_folder'/cellranger_metric_table' \
                                                   --cellranger_long_metrics $experiment_folder'/cellranger_metrics'

/usr/bin/time general_report.R --input $SAMPLES_FILE \
                               --output $report_folder \
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
                               --experiment_name $experiment_name \
                               --results_folder $FULL_RESULTS"/*/preprocessing.R_0000/*" \
                               --resolution $preproc_resolution \
                               --int_sec_cond $int_sec_cond

if [ "$integrative_analysis" == "TRUE" ] ; then
    mkdir -p $FULL_RESULTS/$subset_column
    prior_integration.R --exp_design $exp_design \
                        --output $FULL_RESULTS/$subset_column \
                        --condition $subset_column \
                        --integration_file $integration_file \
                        --experiment_name $experiment_name \
                        --count_path $FULL_RESULTS"/*/cellranger_0000/*" \
                        --suffix "outs/filtered_feature_bc_matrix"
    while IFS= read group; do
    preprocessing.R --input $FULL_RESULTS/$subset_column/$experiment_name"."$group".before.seu.RDS" \
                 --output $FULL_RESULTS/$subset_column \
                 --name $subset_column \
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
                 --report_folder $report_folder \
                 --experiment_name $experiment_name \
                 --resolution $preproc_resolution \
                 --integrate   
    done < $integration_file
fi

