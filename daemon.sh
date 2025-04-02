#! /usr/bin/env bash


framework_dir=`dirname $0`
export CODE_PATH=$(readlink -f $framework_dir )
CONFIG_DAEMON=$1

if [ "$CONFIG_DAEMON" == "" ] ; then
    echo Please specify a config file.
    exit 1
fi

export module=$2 # For setting global vars from config_daemon according to the stage

if [ "$module" == "" ] ; then
    echo Please specify a module to run.
    exit 1
fi

source $CONFIG_DAEMON
mkdir -p $output/report
mkdir -p $output/embeddings
export PATH=$LAB_SCRIPTS:$PATH
export PATH=$CODE_PATH'/aux_parsers:'$PATH
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
    echo Launching sample workflow
    echo $ref_filter > $FULL_RESULTS/ref_filter
    while IFS= read sample; do
        echo Launching $sample
        AF_VARS=`echo "
        \\$sample=$sample,
        \\$exp_design=$exp_design,
        \\$read_path=$read_path,
        \\$FULL_RESULTS=$FULL_RESULTS,
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
        \\$loadRDS=$loadRDS,
        \\$ref_filter=$FULL_RESULTS/ref_filter
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
        \\$FULL_RESULTS=$FULL_RESULTS,
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
        \\$ref_filter=$ref_filter,
        \\$p_adj_cutoff=$p_adj_cutoff,
        \\$verbose=$verbose,
        \\$reduce=$reduce,
        \\$saveRDS=$saveRDS,
        \\$loadRDS=$loadRDS,
        \\$filter_dataset=$filter_dataset
        " | tr -d [:space:]`
        AutoFlow -w $TEMPLATES -V "$AF_VARS" $aux_opt -o $FULL_RESULTS/$sample -v $RESOURCES
        echo Launching pending and failed jobs for $sample
        flow_logger -e $FULL_RESULTS/$sample -w -l -p
    done < $samples_to_process

elif [ "$module" == "2" ] ; then
    echo "Launching stage 2: Sample comparison"
    . ~soft_bio_267/initializes/init_ruby
    cat $FULL_RESULTS/*/metrics > $experiment_folder'/metrics'
    cat $FULL_RESULTS/*/cellranger_metrics > $experiment_folder'/cellranger_metrics'
    create_metric_table.rb $experiment_folder'/metrics' sample $experiment_folder'/metric_table'
    create_metric_table.rb $experiment_folder'/cellranger_metrics' sample $experiment_folder'/cellranger_metric_table'
    compare_samples.R -o $output"/report" -m $experiment_folder'/metric_table' \
                      -l $experiment_folder'/metrics' -e $experiment_name \
                      --cellranger_metrics $experiment_folder'/cellranger_metric_table' \
                      --cellranger_long_metrics $experiment_folder'/cellranger_metrics'
    doublet_files=`find $FULL_RESULTS/*/sc_Hunter.R_0000/ -name doublet_list.txt`
    cat $doublet_files > $experiment_folder/$experiment_name"_doublets.txt"

elif [ "$module" == "3" ] || [ "$module" == "4" ] ; then
    if [ "$module" == "3" ]; then
        echo "Launching stage 3: Single-cell experiment annotation"
        ## if singularity is TRUE, we're launching through singularity image, therefore
        ## sbatch is not available.
        script="$CODE_PATH/aux_sh/annotate_sc.sh"
        if [ "$sketch" == "TRUE" ]; then
            script="$CODE_PATH/singularity/launch_singularity.sh $script"
        fi
    elif [ "$module" == "4" ]; then
        echo "Launching stage 4: DEG analysis"
        . ~soft_bio_267/initializes/init_ruby
        script="$CODE_PATH/aux_sh/sc_Hunter.sh"
        rm -r $TARGETS_FOLDER/*
        mkdir -p $TARGETS_FOLDER
        eval "$generate_targets"
    fi
    if [ "$launch_login" == TRUE ]; then 
        $script
    else
        sbatch $script --cpus-per-task $int_cpu --mem $int_mem
    fi

elif [ "$module" == "5" ] ; then
    # RESULTS PACKAGING
    echo "Creating Single-Cell results pack"
    create_sc_pack.sh
fi
