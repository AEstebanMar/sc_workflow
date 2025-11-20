#! /usr/bin/env bash

source ~soft_bio_267/initializes/init_degenes_hunter

input=$1

clusters_to_enrichment.R -i $data_dir'/clusters_aggregated.txt' -w $enr_cpu -o $output"/report/"$experiment_name"/clust_enr" -f $funsys \
						 -k $gene_keytype -O $organism --size_item $size_item --size_category $size_category --size_edge $size_edge \
						 --node_label $node_label --hilight $highlight --hilight_alpha $highlight_alpha --mode $enr_mode \
						 
