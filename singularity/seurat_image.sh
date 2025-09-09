#! /usr/bin/env bash


source ~soft_bio_267/initializes/init_singularity
script_dir=`dirname "$0"`
export HTMLREPORT_PATH=~aestebanm/dev_R/htmlreportR
export integrate=TRUE
mode=$1
command=$2
CPU=$3

if [ "$mode" != "shell" ] && [ "$mode" != "exec" ] && [ "$mode" != "rebuild" ]; then
	echo "ERROR: Please specify a valid singularity mode. Was \"$mode\""
	echo "Valid modes are:"
	echo "	* shell, launches in interactive mode."
	echo "	* exec, executes a command (passed as second argument in string form)"
	echo "	* rebuild, uncompresses seurat_image.sif, starts interactive mode,
						   then rebuilds the image."
	exit 1
fi

if [ "$mode" == "rebuild" ]; then
	singularity build --sandbox seurat_image_dir/ ./seurat_image.sif
	singularity shell --writable seurat_image_dir/
	singularity build seurat_image.sif seurat_image_dir/
	echo Build done!
	exit 0
fi

if [ "$mode" == "exec" ] && [ "$command" == "" ] ; then
	echo "ERROR: EXEC MODE SPECIFIED BUT NO COMMAND SUPPLIED."
	exit 1
fi

if [ "$mode" != "exec" ] && [ "$command" != "" ]; then
	echo "WARNING: NON-EXEC MODE SPECIFIED BUT COMMAND SUPPLIED. IGNORING THE LATTER."
	unset command
fi

singularity $mode --bind /mnt2:/mnt2 $CODE_PATH/singularity/seurat_image.sif $command --cpus $CPU
