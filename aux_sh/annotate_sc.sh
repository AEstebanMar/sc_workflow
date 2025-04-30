#! /usr/bin/env bash


# STAGE 3 EXPERIMENT ANNOTATION

#SBATCH -J annotate_sc.sh
#SBATCH --cpus-per-task=52
#SBATCH --mem='182gb'
#SBATCH --constraint=sd
#SBATCH --time=3-20:00:00
#SBATCH --error=job.annot.%J.err
#SBATCH --output=job.annot.%J.out

# Setup

if [ -z $SLURM_CPUS_PER_TASK ]; then
    SLURM_CPUS_PER_TASK=1
fi
OMP_NUM_THREADS=1
OPENBLAS_NUM_THREADS=1
MKL_NUM_THREADS=1
echo OMP_NUM_THREADS IS $OMP_NUM_THREADS
echo OPENBLAS_NUM_THREADS IS $OPENBLAS_NUM_THREADS
echo MKL_NUM_THREADS IS $MKL_NUM_THREADS
echo SLURM_CPUS_PER_TASK IS $SLURM_CPUS_PER_TASK

. ~aestebanm/initializes/init_Hunter_dev

hostname
mkdir -p $output"/report"

exp_doublet_file="$experiment_folder/$experiment_name'_doublets.txt'"

annotate_sc.R --output $output \
              --name $experiment_name \
              --doublet_file "$exp_doublet_file" \
              --filter $preproc_filter \
              --mincells $preproc_init_min_cells \
              --minfeats $preproc_init_min_feats \
              --minqcfeats $preproc_qc_min_feats \
              --percentmt $preproc_max_percent_mt \
              --normalmethod $preproc_norm_method \
              --scalefactor $preproc_scale_factor \
              --hvgs $preproc_select_hvgs \
              --ndims $preproc_pca_n_dims \
              --dimheatmapcells $preproc_pca_n_cells \
              --resolution $preproc_resolution \
              --exp_design $exp_design \
              --subset_by "$subset_by" \
              --input $FULL_RESULTS"/*/cellranger_0000/*" \
              --suffix "outs/filtered_feature_bc_matrix" \
              --cluster_annotation "$cluster_annotation" \
              --target_genes $target_genes \
              --cpu $SLURM_CPUS_PER_TASK \
              --imported_counts "$imported_counts" \
              --cell_annotation "$cell_annotation" \
              --SingleR_ref "$refs_path/$SingleR_ref" \
              --ref_version "$ref_version" \
              --ref_label "$ref_label" \
              --ref_de_method "$ref_de_method" \
              --ref_n "$ref_n" \
              --ref_filter "$ref_filter" \
              --p_adj_cutoff $p_adj_cutoff \
              --verbose $verbose \
              --reduce $reduce \
              --samples_to_integrate $samples_to_process \
              --integrate TRUE \
              --int_method "$int_method" \
              --filter_dataset "$filter_dataset" \
              --sketch $sketch \
              --sketch_pct $sketch_pct \
              --sketch_method "$sketch_method"  \
              --force_ncells "$force_ncells" \
              --extra_columns "$extra_columns" \
              --k_weight $k_weight \
              --genome $genome #& process_monitoring.sh R $output/exec_params
if [ ! -s $output/counts/features.tsv.gz ]; then
    echo "Compressing counts data..."
    gzip $output/counts/barcodes.tsv $output/counts/genes.tsv $output/counts/matrix.mtx
    mv $output/counts/genes.tsv.gz $output/counts/features.tsv.gz
fi
