#Author Spealman P. 2022
#
# #Qiime2 Invocation
	# # Version 3 for linux laptop
	## conda activate qiime2-2023.2
#
# # Script Section
	# #
	# #Define directory
	base_dir=/scratch/ps163/Project_Osun/
	mapping_file=${base_dir}/metadata/MappingFile_Osun.tsv
	fastq_dir=${base_dir}/fastq/
	qiime_results=${base_dir}/qiime_results/
	mkdir -p ${qiime_results}
	# #
# # # DADA2 Section
	# ###########
	# # Handles P site - previous site from Mangrove study
	# ###########
	
	qiime tools import \
	--type 'SampleData[SequencesWithQuality]' \
	--input-path /scratch/ps163/Project_Osun/fastq/P_Mangrove/ \
	--input-format CasavaOneEightSingleLanePerSampleDirFmt \
	--output-path /scratch/ps163/Project_Osun/fastq/P/demux-joined.qza
	
	qiime dada2 denoise-single \
	--i-demultiplexed-seqs /scratch/ps163/Project_Osun/fastq/P/demux-joined.qza \
	--p-trim-left 3 \
	--p-trunc-len 0 \
	--o-representative-sequences /scratch/ps163/Project_Osun/fastq/P/pre-otu-rep-seqs-dada2.qza \
	--o-table /scratch/ps163/Project_Osun/fastq/P/pre-otu-table-dada2.qza \
	--o-denoising-stats /scratch/ps163/Project_Osun/fastq/P/stats-dada2.qza
	
	qiime metadata tabulate \
	--m-input-file /scratch/ps163/Project_Osun/fastq/P/stats-dada2.qza \
	--o-visualization /scratch/ps163/Project_Osun/fastq/P/stats-dada2.qzv

	# ##########
	# #Handles V, S sites - From Current Study
	# #/scratch/ps163/Project_Osun/fastq/MSV/casava-18-paired-end-demultiplexed
	# ##########
	
	qiime tools import \
	  --type 'SampleData[PairedEndSequencesWithQuality]' \
	  --input-path /scratch/ps163/Project_Osun/fastq/SV/ \
	  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
	  --output-path /scratch/ps163/Project_Osun/fastq/SV/demux-paired-end.qza
	
	qiime dada2 denoise-paired \
	  --i-demultiplexed-seqs /scratch/ps163/Project_Osun/fastq/SV/demux-paired-end.qza \
	  --p-trim-left-f 13 \
	  --p-trim-left-r 13 \
	  --p-trunc-len-f 150 \
	  --p-trunc-len-r 150 \
	  --o-table /scratch/ps163/Project_Osun/fastq/SV/pre-otu-table-dada2.qza \
	  --o-representative-sequences /scratch/ps163/Project_Osun/fastq/SV/pre-otu-rep-seqs-dada2.qza \
	  --o-denoising-stats /scratch/ps163/Project_Osun/fastq/SV/denoising-stats.qza
	  
	qiime metadata tabulate \
	--m-input-file /scratch/ps163/Project_Osun/fastq/SV/denoising-stats.qza \
	--o-visualization /scratch/ps163/Project_Osun/fastq/SV/denoising-stats.qzv
	  
	###
	# Combine
	###
	#!!!!!!
	qiime feature-table merge \
	--i-tables /scratch/ps163/Project_Osun/fastq/P/pre-otu-table-dada2.qza \
		/scratch/ps163/Project_Osun/fastq/SV/pre-otu-table-dada2.qza \
	--o-merged-table ${qiime_results}/table-dada2.qza \
	
	qiime feature-table merge-seqs \
	--i-data /scratch/ps163/Project_Osun/fastq/P/pre-otu-rep-seqs-dada2.qza \
		/scratch/ps163/Project_Osun/fastq/SV/pre-otu-rep-seqs-dada2.qza \
	--o-merged-data ${qiime_results}/rep-seqs-dada2.qza
	
	# # #v_2:Phylogeny
	qiime phylogeny align-to-tree-mafft-fasttree \
	  --i-sequences ${qiime_results}/rep-seqs-dada2.qza \
	  --o-alignment  ${qiime_results}/aligned-rep-seqs.qza \
	  --o-masked-alignment  ${qiime_results}/masked-aligned-rep-seqs.qza \
	  --o-tree  ${qiime_results}/unrooted-tree.qza \
	  --o-rooted-tree  ${qiime_results}/rooted-tree.qza
	
	# # #v_2:to view trees 
	qiime tools export \
	  --input-path ${qiime_results}/unrooted-tree.qza \
	  --output-path ${qiime_results}/exported-unrooted-tree
	  
	# #v_2:to view trees 
	qiime tools export \
	  --input-path ${qiime_results}/rooted-tree.qza \
	  --output-path ${qiime_results}/exported-rooted-tree
	  
	qiime feature-table summarize \
	  --i-table ${qiime_results}/table-dada2.qza \
	  --o-visualization ${qiime_results}/table-dada2.qzv \
	  --m-sample-metadata-file ${mapping_file}

	qiime feature-table tabulate-seqs \
	  --i-data ${qiime_results}/rep-seqs-dada2.qza \
	  --o-visualization ${qiime_results}/rep-seqs.qzv

	# # ##  
	# # Stop here -
	# # 1. Check min_depth and max_depth
		# # from ${qiime_results}/feature_table.qzv
	# # 2. Using 13000 as minimum:
		# # Retained 104,000 (15.80%) features in 8 (100.00%) samples at the specifed sampling depth.
	# # 3. edit variables below
	# # ##
	# # Alpha - beta section
	# # Alpha raref.
	# # reapeated for 20K, 40K, and 10K for resolution of lower abundant samples
	max_depth=147000
	qiime diversity alpha-rarefaction \
	  --i-table ${qiime_results}/table-dada2.qza \
	  --i-phylogeny ${qiime_results}/rooted-tree.qza \
	  --p-max-depth ${max_depth} \
	  --m-metadata-file ${mapping_file} \
	  --o-visualization ${qiime_results}/alpha-rarefaction_147K.qzv
	  
	max_depth=50000
	qiime diversity alpha-rarefaction \
	  --i-table ${qiime_results}/table-dada2.qza \
	  --i-phylogeny ${qiime_results}/rooted-tree.qza \
	  --p-max-depth ${max_depth} \
	  --m-metadata-file ${mapping_file} \
	  --o-visualization ${qiime_results}/alpha-rarefaction_50K.qzv
	  
	max_depth=13000
	qiime diversity alpha-rarefaction \
	  --i-table ${qiime_results}/table-dada2.qza \
	  --i-phylogeny ${qiime_results}/rooted-tree.qza \
	  --p-max-depth ${max_depth} \
	  --m-metadata-file ${mapping_file} \
	  --o-visualization ${qiime_results}/alpha-rarefaction_13K.qzv
	
	#core diversities
	#double check sampling depth
	min_depth=13000
	core_metrics_results=${qiime_results}/core-metrics-results_5K/
	rm -rf ${core_metrics_results}
		qiime diversity core-metrics-phylogenetic \
		  --i-phylogeny ${qiime_results}/rooted-tree.qza \
		  --i-table ${qiime_results}/table-dada2.qza \
		  --p-sampling-depth ${min_depth} \
		  --m-metadata-file ${mapping_file} \
		  --output-dir ${core_metrics_results}
		  
		#double check sampling depth
		
		qiime diversity alpha-group-significance \
		  --i-alpha-diversity ${core_metrics_results}/observed_features_vector.qza \
		  --m-metadata-file ${mapping_file} \
		  --o-visualization ${core_metrics_results}/observed-water-group-significance.qzv
		
		qiime diversity alpha-group-significance \
		  --i-alpha-diversity ${core_metrics_results}/faith_pd_vector.qza \
		  --m-metadata-file ${mapping_file} \
		  --o-visualization ${core_metrics_results}/faith-pd-group-significance.qzv
		
		qiime diversity alpha-group-significance \
		  --i-alpha-diversity ${core_metrics_results}/shannon_vector.qza \
		  --m-metadata-file ${mapping_file} \
		  --o-visualization ${core_metrics_results}/shannon-group-significance.qzv
		 
		 qiime diversity alpha-group-significance \
		  --i-alpha-diversity ${core_metrics_results}/evenness_vector.qza \
		  --m-metadata-file ${mapping_file} \
		  --o-visualization ${core_metrics_results}/evenness-group-significance.qzv
		
		qiime diversity beta-group-significance \
		  --i-distance-matrix ${core_metrics_results}/weighted_unifrac_distance_matrix.qza \
		  --m-metadata-file ${mapping_file} \
		  --m-metadata-column Site \
		  --o-visualization ${core_metrics_results}/weighted-unifrac-site-significance.qzv \
		  --p-pairwise
		  
		qiime diversity beta-group-significance \
		  --i-distance-matrix  ${core_metrics_results}/bray_curtis_distance_matrix.qza \
		  --m-metadata-file ${mapping_file} \
		  --m-metadata-column Site \
		  --o-visualization ${core_metrics_results}/bray_curtis-site-significance.qzv \
		  --p-pairwise
		  
		qiime diversity beta-group-significance \
		  --i-distance-matrix  ${core_metrics_results}/jaccard_distance_matrix.qza \
		  --m-metadata-file ${mapping_file} \
		  --m-metadata-column Site \
		  --o-visualization ${core_metrics_results}/jaccard_site-significance.qzv \
		  --p-pairwise
		  
		qiime diversity beta-group-significance \
		  --i-distance-matrix  ${core_metrics_results}/jaccard_distance_matrix.qza \
		  --m-metadata-file ${mapping_file} \
		  --m-metadata-column Site \
		  --o-visualization ${core_metrics_results}/jaccard_site-significance.qzv \
		  --p-pairwise
		  
		qiime diversity beta-group-significance \
		  --i-distance-matrix  ${core_metrics_results}/unweighted_unifrac_distance_matrix.qza \
		  --m-metadata-file ${mapping_file} \
		  --m-metadata-column Site \
		  --o-visualization ${core_metrics_results}/unweighted-unifrac-site-significance.qzv \
		  --p-pairwise
# #
# # ##
	#taxonomy
	qiime feature-classifier classify-sklearn \
	  --i-classifier ${base_dir}/metadata/silva-138-99-515-806-nb-classifier.qza \
	  --i-reads ${qiime_results}/rep-seqs-dada2.qza \
	  --o-classification ${qiime_results}/taxonomy.qza
	 
	qiime metadata tabulate \
	  --m-input-file ${qiime_results}/taxonomy.qza \
	  --o-visualization ${qiime_results}/taxonomy.qzv
	  
	qiime taxa barplot \
	  --i-table ${qiime_results}/table-dada2.qza \
	  --i-taxonomy ${qiime_results}/taxonomy.qza \
	  --m-metadata-file ${mapping_file} \
	  --o-visualization ${qiime_results}/taxa-bar-plots.qzv

# ### for genus level 
qiime taxa collapse \
  --i-table ${qiime_results}/table-dada2.qza \
  --i-taxonomy ${qiime_results}/taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table ${qiime_results}/collapsed-genus-table.qza

	qiime taxa barplot \
	  --i-table ${qiime_results}/collapsed-genus-table.qza \
	  --i-taxonomy ${qiime_results}/taxonomy.qza \
	  --m-metadata-file ${mapping_file} \
	  --o-visualization ${qiime_results}/collapsed_genus-taxa-bar-plots.qzv
	
	mkdir extracted-feature-table
	qiime tools extract \
	  --input-path ${qiime_results}/collapsed-genus-table.qza \ \
	  --output-path ${qiime_results}/extracted-feature-table
	 
	qiime tools export \
	--input-path ${qiime_results}/collapsed-genus-table.qza \
	--output-path ${qiime_results}/extracted-feature-table
#
#Particular Alpha Beta tests:
#
	min_depth=5000
	core_metrics_results=${qiime_results}/core-metrics-results_water_special/
	rm -rf ${core_metrics_results}

	
	
	typeis1=water
	typeis2=outflow
	descis=above
#
		qiime feature-table filter-samples \
		  --i-table ${qiime_results}/table-dada2.qza \
		  --m-metadata-file ${mapping_file} \
		  --p-where "([Type]='${typeis1}') OR ([Type]='${typeis2}')" \
		  --o-filtered-table ${qiime_results}/out.qza
		  
		qiime diversity core-metrics-phylogenetic \
		  --i-phylogeny ${qiime_results}/rooted-tree.qza \
		  --i-table ${qiime_results}/out.qza \
		  --p-sampling-depth ${min_depth} \
		  --m-metadata-file ${mapping_file} \
		  --output-dir ${core_metrics_results}
		  
		qiime diversity alpha-group-significance \
		  --i-alpha-diversity ${core_metrics_results}/observed_features_vector.qza \
		  --m-metadata-file ${mapping_file} \
		  --o-visualization ${core_metrics_results}/observed-water-group-significance.qzv
		  
		qiime diversity alpha-group-significance \
		  --i-alpha-diversity ${core_metrics_results}/shannon_vector.qza \
		  --m-metadata-file ${mapping_file} \
		  --o-visualization ${core_metrics_results}/shannon-group-significance.qzv
		  
	min_depth=5000
	core_metrics_results=${qiime_results}/core-metrics-results_sediment_special/
	rm -rf ${core_metrics_results}

	typeis1=sediment
	typeis2=outflow
	descis=above
#
	rm -rf ${core_metrics_results}
		qiime feature-table filter-samples \
		  --i-table ${qiime_results}/table-dada2.qza \
		  --m-metadata-file ${mapping_file} \
		  --p-where "([Type]='${typeis1}') OR ([Type]='${typeis2}')" \
		  --o-filtered-table ${qiime_results}/out.qza
		  
		  
		qiime diversity core-metrics-phylogenetic \
		  --i-phylogeny ${qiime_results}/rooted-tree.qza \
		  --i-table ${qiime_results}/out.qza \
		  --p-sampling-depth ${min_depth} \
		  --m-metadata-file ${mapping_file} \
		  --output-dir ${core_metrics_results}
		  
		# #double check sampling depth
		qiime diversity alpha-group-significance \
		  --i-alpha-diversity ${core_metrics_results}/observed_features_vector.qza \
		  --m-metadata-file ${mapping_file} \
		  --o-visualization ${core_metrics_results}/observed-sediment-group-significance.qzv
		  
		
		qiime diversity alpha-group-significance \
		  --i-alpha-diversity ${core_metrics_results}/shannon_vector.qza \
		  --m-metadata-file ${mapping_file} \
		  --o-visualization ${core_metrics_results}/shannon-group-significance.qzv
		 
		qiime diversity beta-group-significance \
		  --i-distance-matrix ${core_metrics_results}/weighted_unifrac_distance_matrix.qza \
		  --m-metadata-file ${mapping_file} \
		  --m-metadata-column Date \
		  --o-visualization ${core_metrics_results}/weighted-unifrac-site-significance.qzv \
		  --p-pairwise
#
# ANCOM-BC for differential analysis
#
	# #taxonomy
	typeis1=Spring
	typeis2=Valley
		qiime feature-table filter-samples \
		  --i-table ${qiime_results}/table-dada2.qza \
		  --m-metadata-file ${mapping_file} \
		  --p-where "([Site]='${typeis1}') OR ([Site]='${typeis2}')" \
		  --o-filtered-table ${qiime_results}/out.qza
		#
		qiime feature-table summarize\
		  --i-table ${qiime_results}/out.qza \
		  --m-sample-metadata-file ${mapping_file} \
		  --o-visualization ${qiime_results}/out.qzv
		#
		qiime taxa collapse \
		--i-table ${qiime_results}/out.qza \
		--i-taxonomy ${qiime_results}/taxonomy.qza \
		--p-level 6 \
		--o-collapsed-table ${qiime_results}/out_collapsed-genus-table.qza
		#
		qiime composition ancombc \
		--i-table ${qiime_results}/out_collapsed-genus-table.qza \
		--m-metadata-file ${mapping_file} \
		--p-p-adj-method 'fdr' \
		--p-formula Site \
		--o-differentials ${qiime_results}/ancom_bc_${typeis1}_${typeis2}.qza
		
	# #taxonomy
	typeis1=Valley
	typeis2=Mangrove
		qiime feature-table filter-samples \
		  --i-table ${qiime_results}/table-dada2.qza \
		  --m-metadata-file ${mapping_file} \
		  --p-where "([Site]='${typeis1}') OR ([Site]='${typeis2}')" \
		  --o-filtered-table ${qiime_results}/out.qza
		#
		qiime feature-table summarize\
		  --i-table ${qiime_results}/out.qza \
		  --m-sample-metadata-file ${mapping_file} \
		  --o-visualization ${qiime_results}/out.qzv
		#
		qiime taxa collapse \
		--i-table ${qiime_results}/out.qza \
		--i-taxonomy ${qiime_results}/taxonomy.qza \
		--p-level 6 \
		--o-collapsed-table ${qiime_results}/out_collapsed-genus-table.qza
		#
		qiime composition ancombc \
		--i-table ${qiime_results}/out_collapsed-genus-table.qza \
		--m-metadata-file ${mapping_file} \
		--p-p-adj-method 'fdr' \
		--p-formula Site \
		--o-differentials ${qiime_results}/ancom_bc_${typeis1}_${typeis2}.qza
		
		
	# #taxonomy
	typeis1=Spring
	typeis2=Mangrove
		qiime feature-table filter-samples \
		  --i-table ${qiime_results}/table-dada2.qza \
		  --m-metadata-file ${mapping_file} \
		  --p-where "([Site]='${typeis1}') OR ([Site]='${typeis2}')" \
		  --o-filtered-table ${qiime_results}/out.qza
		#
		qiime feature-table summarize\
		  --i-table ${qiime_results}/out.qza \
		  --m-sample-metadata-file ${mapping_file} \
		  --o-visualization ${qiime_results}/out.qzv
		#
		qiime taxa collapse \
		--i-table ${qiime_results}/out.qza \
		--i-taxonomy ${qiime_results}/taxonomy.qza \
		--p-level 6 \
		--o-collapsed-table ${qiime_results}/out_collapsed-genus-table.qza
		#
		qiime composition ancombc \
		--i-table ${qiime_results}/out_collapsed-genus-table.qza \
		--m-metadata-file ${mapping_file} \
		--p-p-adj-method 'fdr' \
		--p-formula Site \
		--o-differentials ${qiime_results}/ancom_bc_${typeis1}_${typeis2}.qza
		