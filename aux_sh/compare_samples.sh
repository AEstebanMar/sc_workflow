#! /usr/bin/env bash


# STAGE 2 SAMPLES COMPARISON

#SBATCH -J compare_samples.sh
#SBATCH --cpus-per-task=12
#SBATCH --mem='1000gb'
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

if [ $imported_counts != "" ]; then

    /usr/bin/time $CODE_PATH/scripts/compare_samples.R -o $report_folder \
                                                       -m $experiment_folder'/metric_table' \
                                                       -l $experiment_folder'/metrics' \
                                                       -e $experiment_name \
                                                       --cellranger_metrics $experiment_folder'/cellranger_metric_table' \
                                                       --cellranger_long_metrics $experiment_folder'/cellranger_metrics'

fi


mkdir -p $FULL_RESULTS/$experiment_name
integration.R --output $FULL_RESULTS/$experiment_name \
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
              --report_folder $report_folder \
              --experiment_name $experiment_name \
              --resolution $preproc_resolution \
              --exp_design $exp_design \
              --int_columns $int_columns \
              --count_path $FULL_RESULTS"/*/cellranger_0000/*" \
              --suffix "outs/filtered_feature_bc_matrix" \
              --samples_to_integrate "$samples_to_integrate" \
              --annotation_dir $annotation_dir \
              --target_genes $exp_data_folder/markers \
              --cpu $SLURM_CPUS_PER_TASK \
              --imported_counts $imported_counts
