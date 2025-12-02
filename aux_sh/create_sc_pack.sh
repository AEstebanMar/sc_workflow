#! /usr/bin/env bash

output_folder=$CODE_PATH/results/$experiment_name

rm -rf $output_folder
mkdir -p $output_folder

DEG_reports=`find $output/report/ -name *DEG_report.html`

cp $output/report/$experiment_name"_All_samples_QC_report.html" $output_folder
cp $output/report/$experiment_name"_annotation_report.html" $output_folder
cp $output/report/$experiment_name"_qc_report.html" $output_folder
cp $DEG_reports $output_folder

for enr_folder in `ls $output/report | grep clust_enr`; do
	mkdir -p $output_folder"/"$enr_folder
	cp -r $output"/report/"$enr_folder/*html $output_folder"/$enr_folder"
done
