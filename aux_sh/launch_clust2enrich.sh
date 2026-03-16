#! /usr/bin/env bash

# STAGE 3 EXPERIMENT ANNOTATION

#SBATCH -J launch_clust2enrich.sh
#SBATCH --constraint=cal
#SBATCH --time=7-00:00:00
#SBATCH --error=job.enr.%J.err
#SBATCH --output=job.enr.%J.out

for TARGET in `ls $output"/DEG"`; do
	# This should be Common_results folder once we introduce compatibility for multiple DEG modes
	DEG_table=$output'/DEG/'$TARGET"/Results*/allgenes*txt"
	func_res_folder=$output"/functional/"$TARGET
	mkdir -p $func_res_folder
	mkdir -p $func_res_folder"/"$enr_subset
	head -n 1 $DEG_table | tr "\t" "\n" > $func_res_folder"/tmp.txt"
	fc_col=`grep "avg_log2FC" $func_res_folder"/tmp.txt"`
	col_vector="cell_type,gene,"$fc_col
	get_columns -i $DEG_table -H -c $col_vector -o $func_res_folder"/gene_fcs"
	cut -f 1,2 $func_res_folder"/gene_fcs" > $func_res_folder"/clusters_file"
	sed -i '1 s/^.*$/cluster\tgeneid\tlog2FC/g' $func_res_folder"/gene_fcs"
	#advanced_options="$advanced_options --gene_attribute_file $func_res_folder/gene_fcs"
	temp_clusters=$func_res_folder"/"$enr_subset"/temp_clusters"
	if [ "$enr_subset" == "global" ]; then
		grep "$enr_subset" $func_res_folder"/clusters_file" > $temp_clusters
		sed -i 1i"Cluster\tgene" $temp_clusters
	elif [ "$enr_subset" == "cell_types" ]; then
		grep -v "global" $func_res_folder"/clusters_file" > $temp_clusters
	fi
	sed -i '1s/cell_type/Cluster/g' $temp_clusters
	input_clusters=$func_res_folder"/"$enr_subset"/"$enr_subset"_clusters.txt"
	aggregate_column_data -i $temp_clusters -x 1 -a 2 -s , | tail -n +2 > $input_clusters
	command="clusters_to_enrichment.R -i $input_clusters -w $enr_cpu \
	-o $output"/report/"$experiment_name"_"$TARGET"_"$enr_subset"_clust_enr" -k $gene_keytype -O $organism \
	--pvalcutoff $enr_pval --qvalcutoff $enr_qval --force $enr_force --clean_parentals $clean_parentals --sim_thr $sim_thr \
	$advanced_resources $advanced_options $advanced_graph_options"
	echo Command called: $command
	/usr/bin/time -o $CODE_PATH"/process_data_single_cell_func" $command
done
wait
