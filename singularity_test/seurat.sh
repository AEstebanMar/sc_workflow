#! /usr/bin/env bash


source ~soft_bio_267/initializes/init_singularity
script_dir=`dirname "$0"`
export singularity=TRUE
export HTMLREPORT_PATH=~aestebanm/dev_R/htmlreportR
mode=$1
command=$2
#singularity build --sandbox seurat/ ./custom_seurat.sif
#singularity shell --writable seurat/
#singularity build custom_seurat.sif seurat/
if [ "$mode" != "shell" ] && [ "$mode" != "exec" ] ; then
	echo "ERROR: Please specify a valid singularity mode. Was $mode"
	exit 1
fi

if [ "$mode" == "exec" ] && [ "$command" == "" ] ; then
	echo "ERROR: EXEC MODE SPECIFIED BUT NO COMMAND SUPPLIED."
	exit 1
fi

if [ "$mode" == "shell" ] && [ "$command" != "" ]; then
	echo "WARNING: SHELL MODE SPECIFIED BUT COMMAND SUPPLIED. IGNORING IT."
	unset command
fi

singularity $mode --bind /mnt2:/mnt2 ./custom_seurat.sif $command

