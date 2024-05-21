#! /usr/bin/env bash

# Sergio Alías, 20230323
# Last modified 20230711

# daemon.sh

# Script for controlling the single-cell Cell Ranger workflow


framework_dir=`dirname $0`
export CODE_PATH=$(readlink -f $framework_dir )
CONFIG_DAEMON=$1
#CONFIG_DAEMON=$CODE_PATH'/config_daemon'
export module=$2 # For setting global vars from config_daemon according to the stage
source $CONFIG_DAEMON
export PATH=$LAB_SCRIPTS:$PATH
export PATH=$CODE_PATH'/scripts:'$PATH
export PATH=$CODE_PATH'/aux_sh:'$PATH
export TEMPLATES=$CODE_PATH'/templates'

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
    while IFS= read sample; do
        echo $sample
            AF_VARS=`echo "
            \\$sample=$sample,
            \\$read_path=$read_path
            " | tr -d [:space:]`
            AutoFlow -w $TEMPLATES"/full_workflow.txt" -V "$AF_VARS" $3 -o $FULL_RESULTS/$sample #$RESOURCES
    done < $SAMPLES_FILE

elif [ "$module" == "2" ] ; then
    # STAGE 2 SAMPLES COMPARISON
    echo "Launching stage 2: Samples comparison"
    if [ $launch_login == TRUE ]; then  
        compare_samples.sh
    else
        sbatch aux_sh/compare_samples.sh
    fi
fi