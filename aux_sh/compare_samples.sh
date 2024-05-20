#! /usr/bin/env bash

# Sergio Alias, 20230530
# Last modified 20230721

# STAGE 2 SAMPLES COMPARISON

#SBATCH -J compare_samples.sh
#SBATCH --cpus-per-task=3
#SBATCH --mem='5gb'
#SBATCH --constraint=cal
#SBATCH --time=0-01:00:00
#SBATCH --error=job.comp.%J.err
#SBATCH --output=job.comp.%J.out

# Setup

. ~soft_bio_267/initializes/init_ruby
. ~soft_bio_267/initializes/init_R
hostname

mkdir -p $report_folder
mkdir -p $PREPROC_RESULTS_FOLDER


cat $FULL_RESULTS/*/metrics > $experiment_folder'/metrics'
cat $FULL_RESULTS/*/cellranger_metrics > $experiment_folder'/cellranger_metrics'
create_metric_table.rb $experiment_folder'/metrics' sample $experiment_folder'/metric_table'
create_metric_table.rb $experiment_folder'/cellranger_metrics' sample $experiment_folder'/cellranger_metric_table'

# Main

/usr/bin/time $CODE_PATH/scripts/compare_samples.R -o $report_folder \
                                                   -m $experiment_folder'/metric_table' \
                                                   -l $experiment_folder'/metrics' \
                                                   -e $experiment_name \
                                                   --cellranger_metrics $experiment_folder'/cellranger_metric_table' \
                                                   --cellranger_long_metrics $experiment_folder'/cellranger_metrics'


if [ "$integrative_analysis" == "TRUE" ] ; then
    export SAMPLES_FILE=$integration_file
    Rscript scripts/prior_integration.R --exp_design $exp_design \
                                        --output $RESULTS_FOLDER \
                                        --condition $subset_column \
                                        --integration_file $integration_file \
                                        --experiment_name $experiment_name \
                                        --count_folder $COUNT_RESULTS_FOLDER
    export SAMPLES_FILE=$integration_file
    preprocessing.R --input run_count)/"$sample"/outs \
                 --output results \
                 --name $sample \
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
                 --report_folder results \
                 --experiment_name $experiment_name \
                 --resolution $preproc_resolution \
                 --integrative_analysis    
fi

# Main

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
                               --results_folder $PREPROC_RESULTS_FOLDER \
                               --resolution $preproc_resolution \
                               --integrative_analysis $integrative_analysis \
                               --int_sec_cond $int_sec_cond