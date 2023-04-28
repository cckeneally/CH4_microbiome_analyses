# Load packages, set up
library(SpiecEasi)
library(devtools)
library(igraph)
library(vegan)
library(Matrix)
library(reshape2)
library(plyr)
library(dplyr)
library(gridExtra)
library(grid)
library(microbiome)
library(phyloseq)
library(qiime2R)
library(tidyverse)
library(microbiomeutilities)

# Data Import
## Set wd
setwd("/home/chriskeneally/Documents/Postgrad/Data/PhoenixOutputs/coorong16sanalysis_2021")

### Import taxonomy, ASVs & metadata
metadata <- read_tsv("sample-metadata.tsv")
ASVs <- read_qza("table-240.qza")
taxonomy <- read_qza("taxonomy-240.qza")
### convert table to tabular split version
taxtable <- taxonomy$data %>%
  as_tibble() %>%
  separate(Taxon, sep=";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))
#make physeq object
phyloseq2021 <- qza_to_phyloseq(
  features = "table-240.qza", 
  tree = "rooted-tree.qza", 
  taxonomy = "taxonomy-240.qza", 
  metadata = "sample-metadata.tsv")
phyloseq2021

#seperate out june samples
phylojune <- subset_samples(phyloseq2021, month=="june")

#summary
microbiome::summarize_phyloseq(phylojune)

# OTU table  
otu_tab <- microbiome::abundances(phylojune)
# check 
otu_tab[1:5,1:5] # for my table show me [1 to 5 rows, 1 to 5 columns]
# Taxonomy table
tax_tab <- phyloseq::tax_table(phylojune)
# check 
tax_tab[1:5,1:5] # for my table show me [1 to 5 otu ids, 1 to 5 first five ranks]

# Clean tax table of uncultured/poorly defined taxa
tax_table(phylojune)[tax_table(phylojune) == "uncultured"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured archaeon"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured euryarchaeote"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured crenarchaeote"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured archaeon 20c-10"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured haloarchaeon"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Halobacteriales archaeon"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured archaeon 2MT16"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Methanosarcinales archaeon"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured archaeon VC2.1 Arc6"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Thermoplasmatales archaeon"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Thermoplasmata archaeon"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "unidentified archaeon"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "unidentified marine bacterioplankton"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured archaeon 20c-39"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured archeon 'KTK 4A'"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured sediment archaeon"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured archaeon 19b-26"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured archaeon 20c-52"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured archaeon TA1f2"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured eukaryote"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured archaeon 19a-29"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured archaeon 19b-39"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured euryarchaeote VAL31-1"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Firmicutes bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured microorganism"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Chloroflexi bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Aminicenantes bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured sediment bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured proteobacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured spirochete"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured cyanobacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "possible genus 03"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured soil bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured candidate division SR1 bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured marine bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured candidate division WS6 bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured anaerobic bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured deep-sea bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Flavobacterium sp."] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Microgenomates group bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured KB1 group bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Dehalococcoides sp."] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Chlorobi bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured delta proteobacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured prokaryote"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured bacterium HF0500_03M05"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Candidatus Gracilibacteria bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Parcubacteria group bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Dehalococcoidia bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Anaerolineae bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Actinomycetales bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uuncultured Clostridia bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Lentisphaerae bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured actinomycete"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured soil bacterium PBS-III-18"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured candidate division BRC1 bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured soil bacterium PBS-III-4"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured actinobacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Latescibacteria bacterium"] <- "Latescibacteria"
tax_table(phylojune)[tax_table(phylojune) == "uncultured bacterium mle1-16"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured soil bacterium PBS-III-30"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Verrucomicrobia bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Acidobacteria bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Omnitrophica bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Kiritimatiellaeota bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Verrucomicrobium sp."] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured planctomycete"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Planctomyces sp."] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Planctomycetales bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Deferribacteres bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Chitinivibrionia bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Caldithrix sp."] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured candidate division TA06 bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Clostridia bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured verrucomicrobium DEV007"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Geobacter sp."] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured candidate division GN04 bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Ignavibacteriales bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Bacteroidetes bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Cytophagales bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Epsilonproteobacteria bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Bacteroidetes/Chlorobi group bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured gamma proteobacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Cytophaga sp."] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured candidate division KSB1 bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Bacteroidales bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Acidobacteriaceae bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Syntrophobacterales bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Actinobacteridae bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured forest soil bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Gemmatimonadetes bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Gloeobacter sp."] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured low G+C Gram-positive bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured alpha proteobacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured archaeon"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Desulfuromonadales bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured compost bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured bacterium zdt-33i5"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Rickettsiales bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Ochrobactrum sp."] <- "Ochrobactrum"
tax_table(phylojune)[tax_table(phylojune) == "uncultured Rhizobiales bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Marinobacter sp."] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured candidate division GN06 bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "Verruc-01"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "Unknown Family"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured bacterium GR-WP33-58"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Polyangiaceae bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Kofleriaceae bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured eubacterium AB16"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Desulfobulbaceae bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Geobacteraceae bacterium"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured Marinobacter sp."] <- ""
tax_table(phylojune)[tax_table(phylojune) == "metagenome"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "marine metagenome"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "hypersaline lake metagenome"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "wastewater metagenome"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "bioreactor metagenome"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "umicrobial mat metagenome"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "anaerobic digester metagenome"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "hydrothermal vent metagenome"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "Ambiguous_taxa"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured archaeon"] <- ""
tax_table(phylojune)[tax_table(phylojune) == "uncultured organism"] <- ""


taxa_names(phylojune) <- paste0("ASV", seq(ntaxa(phylojune)))
#besthit
ps.bh <- format_to_besthit(phylojune)
ps.bh.g <- tax_glom(ps.bh, "Genus")

#subset dep and nondep
ps.dep <- subset_samples(ps.bh, dep == "Depositional")
ps.nodep <- subset_samples(ps.bh, dep == "Other")


# Network/Co-occurence

# load in the following function
'%!in%' <- function(x,y) {
  !('%in%'(x,y))
}

## Infer and test a network - single amplicon
spieceasi.dep.net <- spiec.easi(ps.dep, 
                                method = 'mb',lambda.min.ratio=1e-2, 
                                nlambda=20,
                                icov.select.params = list(rep.num = 50))

# Extract the adjacency matrix from the spiec.easi object. This indicates which pairs of OTUs are adjacent or not in the graph.
spieceasi.matrix <- symBeta(getOptBeta(spieceasi.net), mode='maxabs')
spieceasi.matrix.dsc <- spieceasi.matrix
spieceasi.matrix <- as.matrix(spieceasi.matrix)
# Add the OTU names to the adjacency matrix.
rownames(spieceasi.matrix) <- rownames(tax_table(ps.gen2))
colnames(spieceasi.matrix) <- rownames(tax_table(ps.gen2))
otu.names <- rownames(tax_table(ps.gen2))

# Build a weighted network from the adjacency matrix. The edges in a weighted network represent the strength of association between OTUs.
net <- graph.adjacency(spieceasi.matrix, mode = "undirected",weighted = TRUE, diag = FALSE)
V(net)$name <- otu.names

#Convert the edge weights into distances, where larger weights become shorter distances, and then output a distance-based network
net.dist <- net
max(abs(E(net.dist)$weight))
weights.dist <- 1 - abs(E(net.dist)$weight)
E(net.dist)$weight <- weights.dist

#Convert the weighted network to a separate absolute network
net.abs <- net
E(net.abs)$weight <- abs(E(net.abs)$weight)

#Calculate centrality metrics and create a summary
# alpha centrality
net.alpha <- alpha.centrality(net)
# degree distribution
net.strength <- strength(net.abs)
# betweenness centrality
bet <- betweenness(net.dist,v = V(net.dist))
# make a summary of centrality metrics
summary_cent <- as.data.frame(net.alpha)
colnames(summary_cent) <- ("Alpha_centrality")
rownames(summary_cent) <- rownames(tax_table(ps.gen2))
summary_cent$Weighted_vertex_degree <- net.strength
summary_cent$Betweenness_centrality <- bet
metrics <- summary_cent

#Cluster nodes into modules
wt <- cluster_louvain(net, weights = E(net.dist)$weight)
temp <- V(net)$name
temp <- as.data.frame(temp)
temp$louvain <- membership(wt)
V(net)$louvain <- temp$louvain

#See which nodes have been put into modules with less than or equal to three members and then combine them into a single group that should be consider as not having been assigned a module (group 9).
length(unique(temp$louvain))
summary_modules <- data.frame(table(temp$louvain))
colnames(summary_modules) <- c("louvain", "n")
summary_modules
modules <- as.numeric(summary_modules$louvain[which(summary_modules$n>3)])
x <- max(modules)+1
for (i in c(1:length(temp$temp))) {
  if(temp$louvain[i] %!in% modules){
    temp$louvain[i] <- paste(x)
  }}
modules <- temp
modules$louvain <- as.numeric(modules$louvain)
modules <- modules[order(modules$louvain),]
module.lookup <-
  data.frame("louvain"=unique(modules$louvain),"new_louvain" =
               c(1:length(unique(modules$louvain))))
new <- merge(modules,module.lookup)
modules <- new
modules <- modules[,2:3]
summary_modules <- data.frame(table(modules$new_louvain))
summary_modules
max(modules$new_louvain)

#Test whether the centrality metrics of nodes differ between groups, when considered as a multivariate dataset. Here, we examine whether an OTU’s relative abundance affects its centrality metrics
# to include multiple metrics in the same model they must be z score transformed to all be on the same scale
# z score transformation
metrics.stand <- decostand(metrics, method = "standardize")
# they also cannot be negative so we transform them all to make thempositive
x <- abs(floor(min(metrics.stand)))
metrics.stand.abs <- metrics.stand + x
# Extract the average abundance of each OTU
av.abund <- as.data.frame(rowMeans(otu_table(ps.gen2)))
colnames(av.abund) <- "average_abundance"
# test whether abundance of an OTU significantly influences how important it is in this network
adonis2(metrics.stand.abs ~ av.abund$average_abundance)
#significant

#Differences in individual metrics between groups can also be evaluated.
cor.test(av.abund$average_abundance, metrics$Alpha_centrality, method = "pearson")
#non-sig
cor.test(av.abund$average_abundance, metrics$Weighted_vertex_degree,
         method = "pearson")
#sig
cor.test(av.abund$average_abundance, metrics$Betweenness_centrality,
         method = "pearson")
#sig

#Test whether certain modules have higher centrality than other modules.
# for this we will need to remove the x module category we created

modules.test <- as.data.frame(modules[which(modules$new_louvain !=x),])
colnames(modules.test) <- c("OTU","louvain")
metrics.stand.abs.test <- metrics.stand.abs[which(modules$new_louvain!= x),]
metrics.test <- metrics[which(modules$new_louvain != x),]
# we now test whether modules differ when considering all metrics ina single model
adonis2(metrics.stand.abs.test ~ modules.test$louvain)
#non-sig

#Predict age using only the top ten most central nodes. This step is simply an example of how important nodes can be treated separately
env.all <- data.frame(sample_data(ps.gen2))
otu.table.all <- data.frame(t(otu_table(ps.gen2)))
#repair colnames broken by t()
colnames(otu.table.all) <- rownames(tax_table(ps.gen2))
# env.age <- env.all[-which(is.na(env.all$Replicate)),]
# otu.table.age <- otu.table.all[which(rownames(otu.table.all)%!in%rownames(env.age)),]
adonis2(otu.table.all ~ env.all$Replicate)
#p<0.001
order.alpha <- metrics[order(-metrics$Alpha_centrality),]
top.alpha <- row.names(order.alpha)[1:10]
order.deg <- metrics[order(-metrics$Weighted_vertex_degree),]
top.deg <- row.names(order.deg)[1:10]
order.bet <- metrics[order(-metrics$Betweenness_centrality),]
top.bet <- row.names(order.bet)[1:10]
#subset and remove rows with 0
otu.table.age.alpha <- otu.table.all[,top.alpha]
otu.table.age.deg <- otu.table.all[,top.deg]
otu.table.age.bet <- otu.table.all[,top.bet]
#rows are present with 0 counts, bray cannot be used
adonis2(otu.table.age.alpha ~ env.all$Replicate) #nonsig
adonis2(otu.table.age.deg ~ env.all$Replicate, method = "chisq") #nonsig
adonis2(otu.table.age.bet ~ env.all$Replicate, method = "chisq") #sig

# top.final <- top.deg[which(top.deg %in% top.alpha)]
# otu.table.age.final <- otu.table.age[,top.final]
# top.final
# adonis2(otu.table.age.final ~ env.all$Replicate, method = "chisq")

#Having identified four OTUs that are significantly associated with age we can now export the network to Gephi to view the graph.
# melt the network to prepare for gephi
spieceasi.matrix.m <- reshape2::melt(spieceasi.matrix)
#name cols
colnames(spieceasi.matrix.m) <- c("source","target","weight")
# get the names of all nodes
node.names <-
  unique(c(as.character(unique(spieceasi.matrix.m$source)),as.character
           (unique(spieceasi.matrix.m$target))))
# number them as an alphabetical node list, write to a csv
node.names <- as.data.frame(node.names)
node.names <- as.data.frame(node.names)
node.names$node_number <- c(1:length(node.names$node.names))
node.names$node_number2 <- c(1:length(node.names$node.names))
colnames(node.names) <- c("Taxonomy", "Label", "Id")
row.names(node.names) <- node.names$Taxonomy
row.names(modules) <- modules$temp
modules <- modules[order(modules$temp),]
row.names(node.names) ==row.names(metrics)
row.names(node.names) ==row.names(modules)
node.names.final <- cbind(node.names, metrics,modules)
write.table(node.names.final, "node.names.csv", sep = ",", row.names
            = FALSE)

# create a legend for the network
node.names.label <- data.frame(node.names$Taxonomy, node.names$Label)
colnames(node.names.label) <- c("Taxonomy","Node Label")

library(gridExtra)
g <- tableGrob(node.names.label, rows = NULL)
grid.draw(g)
svg("node_legend.svg", height=, width=4)
grid.draw(g)
dev.off()

#convert node names to numbers, write to a csv
temp <- merge(x = spieceasi.matrix.m, y = node.names, 
              by.x ="source", by.y = "Taxonomy")

#create the edge list
colnames(temp) <-
  c("source","target","weight","remove","source_number")
temp <- temp[,-4]
edge.list <- merge(x = temp, y = node.names, by.x = "target", by.y =
                     "Taxonomy")
colnames(edge.list) <- c("source","target","weight","source.number",
                         "target.number")
edge.list <- edge.list[,c(3,4,6)]
colnames(edge.list) <- c("weight","source","target")
edge.list$Type <- "Undirected"
negative <- ifelse(edge.list$weight<0, "negative", "positive")
edge.list$Negative <- negative
edge.list$weight <- abs(edge.list$weight)
edge.list <- edge.list[which(abs(edge.list$weight)>0),]
write.table(edge.list, "edge_list.csv", sep = ",", row.names = FALSE)
stop()

# EN
# GPr  = transform_sample_counts(phylojune, function(x) x / sum(x) )
# GPfr = filter_taxa(GPr, function(x) mean(x) > 1e-4, TRUE)
# 
# j4 <- make_network(phylojune,
#                    type="taxa", distance="jaccard",
#                    max.dist = 0.4, keep.isolates=FALSE)
# 
# plot_network(j4, methano_june, color = "Replicate", line_weight=0.6) +
#   scale_color_npg()

# co_occurrence_network(methano_abund, 
#                       treatment = 'Replicate', 
#                       subset = c('East Deep','West Deep'),
#                       co_occurrence_table = NULL, 
#                       classification = 'Phylum')

#Inference of Microbial Ecological Networks

#The input for SPIEC-EASI is a counts table. The normalization and tranformation is done by the function.
#This step is heavy on computational memory and slow. Noise filtered OTU-OTU level covariance would be ideal.

#reduce ASVs to test
# ps1.stool.otu <- prune_taxa(taxa_sums(methano_june) > 0, methano_june)
# 
# # Add taxonomic classification to OTU ID
# ps1.stool.otu.f <- microbiomeutilities::format_to_besthit(methano_june)
# 
# ## Warning: replacing previous import 'ggplot2::alpha' by 'microbiome::alpha' when
# ## loading 'microbiomeutilities'
# 
# head(tax_table(ps1.stool.otu))
# 
# 
# #Check the difference in two phyloseq objects.
# 
# head(tax_table(ps1.stool.otu.f))

# Prepare data for SpiecEasi

#The calculation of SpiecEasi are time consuming. For this tutorial, we will have the necessary input files for SpiecEasi.

se.mgen.besthit <- spiec.easi(ps.mgen.besthit, method='mb', lambda.min.ratio=1e-2,
                              nlambda=99)
mgen.ig2.mb <- adj2igraph(getRefit(se.mgen.besthit), 
                          vertex.attr=list(name=taxa_names(ps.mgen.besthit)))
netplot <- plot_network(mgen.ig2.mb, ps.mgen.besthit, type='taxa', color="Genus")

se.ps.gen.beshit <- spiec.easi(ps.gen.besthit, method='mb', lambda.min.ratio=1e-2,
                               nlambda=20, pulsar.params=list(rep.num=50))
ig2.mb <- adj2igraph(getRefit(se.mb.amgut2),  vertex.attr=list(name=taxa_names(ps.gen.besthit)))
netplot <- plot_network(ig2.mb, ps.gen.besthit, type='taxa', color="Phylum")

plot(adj2igraph(getRefit(se.mb.amgut2)), vertex.size=9)

# otu.c <- t(otu_table(ps1.stool.otu.f)@.Data) #extract the otu table from phyloseq object
# 
# tax.c <- as.data.frame(tax_table(ps1.stool.otu.f)@.Data)#extract the taxonomy information
# 
# head(tax.c)

# use this only for first attempt to run it on server to save time
#saveRDS(otu.c, "input_data/stool.otu.c.rds")
#saveRDS(tax.c, "input_data/stool.tax.c.rds")

# SPIEC-EASI network reconstruction

# In practice, use 99+ subsamples (rep.num)
set.seed(3)
net.c <- spiec.easi(otu.c, method='mb', icov.select.params=list(rep.num=99)) 
# reps have to increases for real data

# saveRDS(net.c, "input_data/net.c.rds")

#please use more numebr of rep.num (99 or 999) the paraemters 

## Create graph object and get edge values  


# the PC has low processing power, you can read the otuput created by us present in the input_data folder.
#source("scripts/symBeta.R") # load custom function to get weights.
#net.c <- readRDS("input_data/stool.net.rds")
class(net.c)

## [1] "select"

n.c <- symBeta(getOptBeta(net.c))

#Add names to IDs
#We also add abundance values to vertex (nodes).

colnames(n.c) <- rownames(n.c) <- colnames(otu.c)

vsize <- log2(apply(otu.c, 2, mean)) # add log abundance as properties of vertex/nodes.

# 9.2.1 Prepare data for plotting
library(igraph)

stool.ig <- graph.adjacency(n.c, mode='undirected', add.rownames = TRUE, weighted = TRUE)
stool.ig # we can see all the attributes and weights

## IGRAPH a8db9f0 UNW- 679 2454 -- 
## + attr: name (v/c), TRUE (v/c), weight (e/n)
## + edges from a8db9f0 (vertex names):
##  [1] OTU-9410491526:Bacteroides--OTU-9410491516:Bacteroides   
##  [2] OTU-9410491526:Bacteroides--OTU-9410491518:Bacteroides   
##  [3] OTU-9410491526:Bacteroides--OTU-941049327:Bacteroides    
##  [4] OTU-9410491526:Bacteroides--OTU-941049949:Bacteroides    
##  [5] OTU-9410491526:Bacteroides--OTU-9410491514:Bacteroides   
##  [6] OTU-9410491526:Bacteroides--OTU-9410491513:Bacteroides   
##  [7] OTU-9410491526:Bacteroides--OTU-9410491574:Parasutterella
##  [8] OTU-9410491516:Bacteroides--OTU-9410491522:Bacteroides   
## + ... omitted several edges

#plot(stool.ig)

#set the layout option

# check what is it?
?layout_with_fr

coords.fdr <- layout_with_fr(stool.ig)

#9.2.2 igraph network

E(stool.ig)[weight > 0]$color<-"steelblue" #now color the edges based on their values positive is steelblue
E(stool.ig)[weight < 0]$color<-"orange"  #now color the edges based on their values

plot(stool.ig, layout=coords.fdr, vertex.size = 2, vertex.label.cex = 0.5)

The visualisation can be enhanced using ggnet R package.

stool.net <- asNetwork(stool.ig)
network::set.edge.attribute(stool.net, "color", ifelse(stool.net %e% "weight" > 0, "steelblue", "orange"))

Start adding taxonomic information.

colnames(tax_table(ps1.stool.otu.f))

## [1] "Domain"   "Phylum"   "Class"    "Order"    "Family"   "Genus"    "best_hit"

phyla <- map_levels(colnames(otu.c), from = "best_hit", to = "Phylum", tax_table(ps1.stool.otu.f))
stool.net %v% "Phylum" <- phyla
stool.net %v% "nodesize" <- vsize

9.2.3 Network plot

mycolors <- scale_color_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928"))

p <- ggnet2(stool.net, node.color = "Phylum", 
            label = TRUE, node.size = "nodesize", 
            label.size = 2, edge.color = "color") + guides(color=guide_legend(title="Phylum"), size = FALSE) + mycolors

p 

This is difficult to interpret. One way is to remove nodes that are connected to few other nodes. We can use degree as a network statisitic.

stl.mb <- degree.distribution(stool.ig)
plot(0:(length(stl.mb)-1), stl.mb, ylim=c(0,.35), type='b', 
     ylab="Frequency", xlab="Degree", main="Degree Distributions")

# we will look at only taxa connect more than 10 others
p <- ggnet2(stool.net, node.color = "Phylum", 
            label = TRUE, 
            label.size = 3, edge.color = "color",
            size = "degree", size.min = 10) + guides(color=guide_legend(title="Phylum"), size = FALSE) + mycolors

## size.min removed 521 nodes out of 679

## Scale for 'colour' is already present. Adding another scale for 'colour',
## which will replace the existing scale.

p 

9.3 Network properties

Check for the number of positive and negative edges.

betaMat=as.matrix(symBeta(getOptBeta(net.c)))

# We divide by two since an edge is represented by two entries in the matrix.
positive=length(betaMat[betaMat>0])/2 

negative=length(betaMat[betaMat<0])/2 

total=length(betaMat[betaMat!=0])/2 

9.3.1 Modularity in networks

net.c

## Model: Meinshausen & Buhlmann Graph Estimation (mb)
## selection criterion: stars 
## Graph dimension: 679 
## sparsity level 0.01066118

mod.net <- net.c$refit

colnames(mod.net) <- rownames(mod.net) <- colnames(otu.c)#you can remove this 

vsize <- log2(apply(otu.c, 2, mean))# value we may or may not use as vertex.attribute

stool.ig.mod <- graph.adjacency(mod.net, mode='undirected', add.rownames = TRUE)
plot(stool.ig.mod) # we can see all the attributes and weights

stool.net.mod <- asNetwork(stool.ig.mod)

Set vertex attributes. We can color by phyla and set the size of nodes based on log2 abundance.

phyla <- map_levels(colnames(otu.c), from = "best_hit", to = "Phylum", tax_table(ps1.stool.otu.f))
stool.net.mod %v% "Phylum" <- phyla
stool.net.mod %v% "nodesize" <- vsize

9.3.2 Network plot

mycolors <- scale_color_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928"))

# check the colorpicker in the addins option in RStudio to interactively select color options.  

p <- ggnet2(stool.net.mod, node.color = "Phylum", 
            label = TRUE, node.size = 2, 
            label.size = 2) + guides(color=guide_legend(title="Phylum"), size = FALSE) + mycolors

p 

Identify modularity in networks.

modules =cluster_fast_greedy(stool.ig.mod)

print(modules)

## IGRAPH clustering fast greedy, groups: 12, mod: 0.47
## + groups:
##   $`1`
##    [1] "OTU-9410492646:Bacteroides"          
##    [2] "OTU-9410492641:Bacteroides"          
##    [3] "OTU-9410492645:Bacteroides"          
##    [4] "OTU-9410491981:Parabacteroides"      
##    [5] "OTU-9410491974:Bacteroides"          
##    [6] "OTU-9410491978:Bacteroides"          
##    [7] "OTU-9410491977:Parabacteroides"      
##    [8] "OTU-9410491976:Bacteroides"          
##    [9] "OTU-9410492922:Bacteroides"          
##   + ... omitted several groups/vertices

modularity(modules)

## [1] 0.4725467

V(stool.ig.mod)$color=modules$membership

plot(stool.ig.mod, col = modules, vertex.size = 4, vertex.label = NA)

stool.net.mod %v% "membership" <- modules$membership

p <- ggnet2(stool.net.mod, node.color = "membership", 
            label = TRUE, node.size = "nodesize", 
            label.size = 2) + guides(color=guide_legend(title="membership"), size = FALSE) + mycolors

## Scale for 'colour' is already present. Adding another scale for 'colour',
## which will replace the existing scale.

#Check which OTUs are part of different modules.

modulesOneIndices=which(modules$membership==1)
modulesOneOtus=modules$names[modulesOneIndices]
modulesTwoIndices=which(modules$membership==2)
modulesTwoOtus=modules$names[modulesTwoIndices]

modulesThreeIndices=which(modules$membership==3)
modulesThreeOtus=modules$names[modulesThreeIndices]
modulesFourIndices=which(modules$membership==4)
modulesFourOtus=modules$names[modulesFourIndices]

modulesFiveIndices=which(modules$membership==5)
modulesFiveOtus=modules$names[modulesFiveIndices]
modulesSixIndices=which(modules$membership==6)
modulesSixOtus=modules$names[modulesSixIndices]

print(modulesOneOtus)


#sparcc networks
sparcc.ps <- sparcc(as.matrix(otu_table(ps.bh)), iter = 99)

#Pseudopvals
#psep.ps <- sparccboot(as.matrix(otu_table(ps.bh)), R = 100)
#pval.ps <- pval.sparccboot()

## Define arbitrary threshold for SparCC correlation matrix for the graph
sparcc.dep.bh.graph <- abs(sparcc.ps$Cor) >= 0.1
diag(sparcc.dep.bh.graph) <- 0
library(Matrix)
sparcc.dep.bh.graph <- Matrix(sparcc.dep.bh.graph, sparse=TRUE)
## Create igraph objects
ig.dep.bh <- adj2igraph(sparcc.dep.bh.graph)

library(RCy3)

createNetworkFromIgraph(ig.dep.bh,"myIgraph")

# Visualise
library(igraph)
## set size of vertex proportional to clr-mean
vsize    <- rowMeans(clr(otu_table(ps.bh), 1))+6
am.coord <- layout.fruchterman.reingold(ig.dep.bh)

par(mfrow=c(1,3))
plot(ig.dep.bh, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="sparcc")