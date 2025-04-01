#! /usr/bin/env bash


source ~soft_bio_267/initializes/init_singularity
script_dir=`dirname "$0"`
export HTMLREPORT_PATH=~aestebanm/dev_R/htmlreportR
mode=$1
command=$2

if [ "$mode" == "rebuild" ]; then
	singularity build --sandbox seurat_image_dir/ ./seurat_image.sif
	singularity shell --writable seurat_image_dir/
	singularity build seurat_image.sif seurat_image_dir/
	echo Build done!
	exit 0
fi
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

singularity $mode --bind /mnt2:/mnt2 $CODE_PATH/singularity/seurat_image.sif $command
