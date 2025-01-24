#! /usr/bin/env bash

source ~soft_bio_267/initializes/init_degenes_hunter
output_folder=$CODE_PATH/results/$experiment_name

rm -rf $output_folder
mkdir -p $output_folder
mkdir -p $output_folder

cp $output/report/$experiment_name"_All_samples_QC_report.html" $output_folder
cp $output/report/$experiment_name"_integration_report.html" $output_folder
cp $output/report/$experiment_name"_qc_report.html" $output_folder
cp $output/$experiment_name".final_results.rds" $output_folder
