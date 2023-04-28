## PICRUSt pre-processing
## Starting with decontaminated physeq file
mean(sample_sums(phylojune)) #51912.21 ASVs
phylojune #18432 Taxa


#Filters ASVS with <5 reads, and not present in at least 2 samples
filter <- phyloseq::genefilter_sample(phylojune, filterfun_sample(function(x) x >= 5), A = 2)
pi_phylojune <- phyloseq::prune_taxa(filter, phylojune)

mean(sample_sums(pi_phylojune)) #49609.76 ASVs
pi_phylojune #4467 taxa

library(biomformat)
biom_pi <- make_biom(otu_table(pi_phylojune))
output_path <- "~/Documents/Postgrad/Data/PhoenixOutputs/coorong16sanalysis_2021/PICRUST2/biom_pi.biom"  # Set the output file path
write_biom(biom_pi, output_path)
