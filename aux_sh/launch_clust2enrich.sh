#! /usr/bin/env bash


for TARGET in `ls $output"/DEG"`; do
	# This should be Common_results folder once we introduce compatibility for multiple DEG modes
	DEG_table=$output'/DEG/'$TARGET"/Results*/allgenes*txt"
	func_res_folder=$output"/functional/"$TARGET
	mkdir -p $func_res_folder
	mkdir -p $func_res_folder"/"$enr_subset
	get_columns -i $DEG_table -H -c cell_type,gene -o $func_res_folder"/clusters_file"
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
	/usr/bin/time -o $CODE_PATH"/process_data_single_cell_func" clusters_to_enrichment.R -i $input_clusters -w $enr_cpu \
	-o $output"/report/"$experiment_name"_"$TARGET"_"$enr_subset"_clust_enr" -f $funsys -k $gene_keytype -O $organism --size_item $size_item \
	--size_category $size_category --size_edge $size_edge --node_label $node_label --hilight $highlight --hilight_alpha $highlight_alpha \
	--mode $enr_mode -F
done
wait
