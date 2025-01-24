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
mkdir -p $output/report
export PATH=$LAB_SCRIPTS:$PATH
export PATH=$CODE_PATH'/scripts:'$PATH
export PATH=$CODE_PATH'/aux_sh:'$PATH
export TEMPLATE_PATH=$CODE_PATH'/templates'

aux_opt=$3

if [ "$imported_counts" != "" ]; then
    TEMPLATES=$TEMPLATE_PATH/divide_counts.af
else
    TEMPLATES=$TEMPLATE_PATH/count_sc.af
fi

TEMPLATES="$TEMPLATES,$TEMPLATE_PATH/sc_sample_analysis.af"

. ~soft_bio_267/initializes/init_autoflow

## STAGE EXECUTION
#######################################################################

if [ "$module" == "0" ] ; then
    # STAGE 0: REFERENCE PREPARATION
    echo "Generating reference from specified database"
    if [ -d "$ref_origin" ]; then
        sbatch $CODE_PATH/aux_sh/get_SingleR_ref.sh $aux_opt
    else
        get_SingleR_ref.sh $aux_opt
    fi
fi

if [ "$module" == "1" ] ; then
    mkdir -p $FULL_RESULTS
    echo Launching workflow
    while IFS= read sample; do
        echo Launching $sample
        AF_VARS=`echo "
        \\$sample=$sample,
        \\$exp_design=$exp_design,
        \\$read_path=$read_path,
        \\$aux_sh_dir=$CODE_PATH/aux_sh,
        \\$script_dir=$CODE_PATH/scripts,
        \\$report_folder=$output/report,
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
        \\$target_genes=$target_genes,
        \\$imported_counts=$imported_counts,
        \\$output=$output,
        \\$cell_annotation=$cell_annotation,
        \\$SingleR_ref=$refs_path/$SingleR_ref,
        \\$ref_version=$ref_version,
        \\$ref_label=$ref_label,
        \\$ref_de_method=$ref_de_method,
        \\$ref_n=$ref_n,
        \\$p_adj_cutoff=$p_adj_cutoff,
        \\$verbose=$verbose,
        \\$reduce=$reduce,
        \\$saveRDS=$saveRDS,
        \\$loadRDS=$loadRDS
        " | tr -d [:space:]`
        AutoFlow -w $TEMPLATES -V "$AF_VARS" $aux_opt -o $FULL_RESULTS/$sample $RESOURCES
    done < $samples_to_process

elif [ "$module" == "1b" ] ; then
    echo Checking workflow execution
    while IFS= read sample; do
        echo Sample $sample
        flow_logger -e $FULL_RESULTS/$sample -w -r all
    done < $samples_to_process

elif [ "$module" == "1c" ] ; then
    echo Regenerating code
    while IFS= read sample; do
        AF_VARS=`echo "
        \\$sample=$sample,
        \\$exp_design=$exp_design,
        \\$read_path=$read_path,
        \\$aux_sh_dir=$CODE_PATH/aux_sh,
        \\$script_dir=$CODE_PATH/scripts,
        \\$report_folder=$output/report,
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
        \\$target_genes=$target_genes,
        \\$imported_counts=$imported_counts,
        \\$output=$output,
        \\$cell_annotation=$cell_annotation,
        \\$SingleR_ref=$refs_path/$SingleR_ref,
        \\$ref_version=$ref_version,
        \\$ref_label=$ref_label,
        \\$ref_de_method=$ref_de_method,
        \\$ref_n=$ref_n,
        \\$p_adj_cutoff=$p_adj_cutoff,
        \\$verbose=$verbose,
        \\$reduce=$reduce,
        \\$saveRDS=$saveRDS,
        \\$loadRDS=$loadRDS
        " | tr -d [:space:]`
        AutoFlow -w $TEMPLATES -V "$AF_VARS" $3 -o $FULL_RESULTS/$sample -v $RESOURCES
        echo Launching pending and failed jobs for $sample
        flow_logger -e $FULL_RESULTS/$sample -w -l -p
    done < $samples_to_process

elif [ "$module" == "2" ] ; then
    # STAGE 2 SEURAT ANALYSIS
    echo "Launching stage 2: Samples comparison"
    if [ $launch_login == TRUE ]; then  
        compare_samples.sh
    else
        sbatch $CODE_PATH/aux_sh/compare_samples.sh --cpus-per-task $CPU --mem $mem
    fi
    
elif [ "$module" == "3" ] ; then
    # RESULTS PACKAGING
    echo "Creating Single-Cell results pack"
    create_sc_pack.sh
fi
