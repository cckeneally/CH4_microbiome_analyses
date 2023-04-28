#GHUB .ipynb

import qiime2
from qiime2 import Artifact
from qiime2.plugins.feature_classifier.methods import classify_sklearn,\
                                                      extract_reads, \
                                                      fit_classifier_naive_bayes

# Training feature classifiers with q2-feature-classifier
# https://docs.qiime2.org/2019.1/tutorials/feature-classifier/
silva_132 = Artifact.import_data('FeatureData[Sequence]',
                                 'rep_set/rep_set_16S_only/99/silva_132_99_16S.fna')

silva_132_taxonomy = Artifact.import_data('FeatureData[Taxonomy]',
                                          'taxonomy/16S_only/99/majority_taxonomy_7_levels.txt',
                                           view_type = 'HeaderlessTSVTaxonomyFormat')
   


# extract reference reads
# V3-V4 region: 341f: CCTACGGGNGGCWGCAG; 806r: GACTACHVGGGTATCTAATCC
qiime feature-classifier extract-reads \
  --i-sequences S132_99_sequences.qza \
  --p-f-primer GTGCCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --p-trunc-len 120 \
  --p-min-length 100 \
  --p-max-length 400 \
  --o-reads reference_sequences.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads reference_sequences.qza \
  --i-reference-taxonomy reference_taxonomy.qza \
  --o-classifier classifier.qza