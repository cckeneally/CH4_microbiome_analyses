$ conda activate qiime2-2021.8

#Paired end fastq manifest (Phred 33V2) created to import paired ends into single samples

$qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path se-33-manifest.txt \
--output-path paired-end-demux.qza \
--input-format PairedEndFastqManifestPhred33V2

#Imported se-33-manifest.txt as PairedEndFastqManifestPhred33V2 to paired-end-demux.qza

$qiime demux summarize \
--i-data paired-end-demux.qza \
--o-visualization paired-end-demux.qzv

# Saved Visualization to: paired-end-demux.qzv

#full length with 13bp trim L/R to ensure paired end overlaps (V3-V4 region ~450bp + overlap)

$qiime dada2 denoise-paired \
--i-demultiplexed-seqs paired-end-demux.qza \
--p-trim-left-f 13 \
--p-trim-left-r 13 \
--p-trunc-len-f 300 \
--p-trunc-len-r 300 \
--o-table table.qza \
--o-representative-sequences rep-seqs.qza \
--o-denoising-stats denoising-stats.qza
 
$qiime feature-table summarize \
--i-table table.qza \
--o-visualization table.qzv \
--m-sample-metadata-file sample-metadata.tsv

#table.qzv visualisation informs sampling depth,
#value=10600 chosen to ensure adequate diversity capture, but loss of 2 samples (T1S1_2, T1S2_2 <10600 reads)
#Retained 487,600 (44.95%) features in 46 (95.83%) samples at the specifed sampling depth.

#Alphararefaction helps to visualise alpha diversity @ passed sampling depth, diversity plateus at 10600 reads.

$qiime diversity alpha-rarefaction \
--i-table table.qza \
--i-phylogeny rooted-tree.qza \
--p-max-depth 4000 \
--m-metadata-file sample-metadata.tsv \
--o-visualization alpha-rarefaction.qzv

$qiime feature-table tabulate-seqs \
--i-data rep-seqs.qza \
--o-visualization rep-seqs.qzv

$qiime metadata tabulate \
--m-input-file denoising-stats.qza \
--o-visualization denoising-stats.qzv

#Creating tree for phylogenetic diversity analyses
$qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences rep-seqs.qza \
--o-alignment aligned-rep-seqs.qza \
--o-masked-alignment masked-aligned-rep-seqs.qza \
--o-tree unrooted-tree.qza \
--o-rooted-tree rooted-tree.qza

#see line 36 comment for --p-sampling-depth rationale

$qiime diversity core-metrics-phylogenetic \
--i-phylogeny rooted-tree.qza \
--i-table table.qza \
--p-sampling-depth 10600 \ 
--m-metadata-file sample-metadata.tsv \
--output-dir core-metrics-results
 
#Weighted taxonomic classifier used cite:
#Michael S Robeson II, Devon R O’Rourke, Benjamin D Kaehler, Michal Ziemski, Matthew R Dillon, Jeffrey T Foster, Nicholas A Bokulich. RESCRIPt: Reproducible sequence taxonomy reference database management for the masses. bioRxiv 2020.10.05.326504; doi: https://doi.org/10.1101/2020.10.05.326504
#Bokulich, N.A., Kaehler, B.D., Rideout, J.R. et al. Optimizing taxonomic classification of marker-gene amplicon sequences with QIIME 2’s q2-feature-classifier plugin. Microbiome 6, 90 (2018). https://doi.org/10.1186/s40168-018-0470-z
#Kaehler, B.D., Bokulich, N.A., McDonald, D. et al. Species abundance information improves sequence taxonomy classification accuracy. Nature Communications 10, 4643 (2019). https://doi.org/10.1038/s41467-019-12669-6

#Tax analysis with classifier based on weighted Silva 139 99% fullseqs 
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-nb-weighted-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv


