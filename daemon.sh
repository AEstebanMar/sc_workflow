#! /usr/bin/env bash

# Sergio Alías, 20230323
# Last modified 20230711

# daemon.sh

# Script for controlling the single-cell Cell Ranger workflow


framework_dir=`dirname $0`
export CODE_PATH=$(readlink -f $framework_dir )
CONFIG_DAEMON=$1
export module=$2 # For setting global vars from config_daemon according to the stage
source $CONFIG_DAEMON
export PATH=$LAB_SCRIPTS:$PATH
export PATH=$CODE_PATH'/scripts:'$PATH
export PATH=$CODE_PATH'/aux_sh:'$PATH
export TEMPLATE_PATH=$CODE_PATH'/templates'

if [ "$imported_counts" != "" ]; then
    TEMPLATES=$TEMPLATE_PATH/divide_counts.af
else
    TEMPLATES=$TEMPLATE_PATH/count_sc.af
fi

TEMPLATES="$TEMPLATES,$TEMPLATE_PATH/preprocessing.af"

. ~soft_bio_267/initializes/init_autoflow

## STAGE EXECUTION
#######################################################################

# if [ "$module" == "0" ] ; then
#     # STAGE 0 CONVERTING BCL FILES INTO FASTQ
#     echo "Launching stage 0: Converting BCL files into FASTQ"
#     if [ $launch_login == TRUE ]; then  
#         cellranger_mkfastq.sh
#     else
#         sbatch aux_sh/cellranger_mkfastq.sh
#     fi
# el
if [ "$module" == "1" ] ; then
    mkdir -p $FULL_RESULTS
    echo Launching workflow
    while IFS= read sample; do
        echo Launching $sample
        AF_VARS=`echo "
        \\$sample=$sample,
        \\$read_path=$read_path,
        \\$aux_sh_dir=$CODE_PATH/aux_sh,
        \\$script_dir=$CODE_PATH/scripts,
        \\$preproc_filter=$preproc_filter,
        \\$preproc_init_min_cells=$preproc_init_min_cells,
        \\$preproc_init_min_feats=$preproc_init_min_feats,
        \\$preproc_qc_min_feats=$preproc_qc_min_feats,
        \\$preproc_max_percent_mt=$preproc_max_percent_mt,
        \\$preproc_norm_method=$preproc_norm_method,
        \\$preproc_scale_factor=$preproc_scale_factor,
        \\$preproc_select_hvgs=$preproc_select_hvgs,
        \\$preproc_pca_n_dims=$preproc_pca_n_dims,
        \\$preproc_pca_n_cells=$preproc_pca_n_cells,
        \\$experiment_name=$experiment_name,
        \\$preproc_resolution=$preproc_resolution,
        \\$imported_counts=$imported_counts
        " | tr -d [:space:]`
        AutoFlow -w $TEMPLATES -V "$AF_VARS" $3 -o $FULL_RESULTS/$sample
    done < $SAMPLES_FILE

elif [ "$module" == "1b" ] ; then
    echo Checking workflow execution
    while IFS= read sample; do
        flow_logger -e $FULL_RESULTS/$sample -w -r all
    done < $SAMPLES_FILE

elif [ "$module" == "1c" ] ; then
    echo Regenerating code
    while IFS= read sample; do
        AF_VARS=`echo "
        \\$sample=$sample,
        \\$read_path=$read_path,
        \\$aux_sh_dir=$CODE_PATH/aux_sh,
        \\$script_dir=$CODE_PATH/scripts,
        \\$preproc_filter=$preproc_filter,
        \\$preproc_init_min_cells=$preproc_init_min_cells,
        \\$preproc_init_min_feats=$preproc_init_min_feats,
        \\$preproc_qc_min_feats=$preproc_qc_min_feats,
        \\$preproc_max_percent_mt=$preproc_max_percent_mt,
        \\$preproc_norm_method=$preproc_norm_method,
        \\$preproc_scale_factor=$preproc_scale_factor,
        \\$preproc_select_hvgs=$preproc_select_hvgs,
        \\$preproc_pca_n_dims=$preproc_pca_n_dims,
        \\$preproc_pca_n_cells=$preproc_pca_n_cells,
        \\$experiment_name=$experiment_name,
        \\$preproc_resolution=$preproc_resolution,
        \\$imported_counts=$imported_counts
        " | tr -d [:space:]`
        AutoFlow -w $TEMPLATES -V "$AF_VARS" $3 -o $FULL_RESULTS/$sample -v
        echo Launching pending and failed jobs for $sample
        flow_logger -e $FULL_RESULTS/$sample -w -l -p
    done < $SAMPLES_FILE

elif [ "$module" == "2" ] ; then
    # STAGE 2 SEURAT ANALYSIS
    echo "Launching stage 2: Samples comparison"
    if [ $launch_login == TRUE ]; then  
        compare_samples.sh
    else
        sbatch aux_sh/compare_samples.sh
    fi
fi