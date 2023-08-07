Coorong 16S Analyses (2021)
================
Christopher Keneally
02/11/2021

\#opts_knit\$set(root.dir = “home/chriskeneally/Desktop/RWD”)

\#Setup \## Load Libraries and other setup

## Set theme

# Set plotting theme

# Analysis (June Sampling)

# Data Import

### Import taxonomy, ASVs & metadata

\#metadata boxplots

\#2

## Metadata Corr

\#Mult regression to find drivers of CH4

# Pairwaise tests for figure 2

## Time Series vars

## Create Phyloseq Object

\#contams

## Data Summarisation

### Sample sequencing depths

\#Sample read count stats

\#Table 1

## Clean taxa & subset into seperate Phyloseq objects

\#June Sample read count stats

# New barplots - Arch

\#New barplots - Bact

## Kingdom level diffabund testing

## Subsets methanogen & Tax analyses

## Subsets Methanotroph & Taxa level analyses

## Physeq tut

\#`{r message=FALSE} #library("ape") #random_tree <- rtree(ntaxa(phylojune), rooted=TRUE, tiplabels=taxa_names(phylojune)) #`

## Taxa bar plots (June - King Level agglom)

### Set up data to visualise relative abundance and group phyla

### Plot (King - June)

Taxa Stacked

### Before going further: Rarefy to even read count depth

### Multiplot: ENV variables ~ Richness (Chao1)

\#first need to check normal data etc

### Uncononstrained ordinations (PCoA and NMDS)

#### Mirlyn - analysis of diversity with multiple iterations of rarefication

#### NMDS Plot - All microbes

### Beta diversity: Homogeneity of group dispersions and compositional dissimilarity/similarity

betadisper - calculates avg distance of group members to group centroid
and ANOVA then tests agains a null hypothesis of: No difference in
dispersion btw groups. ADONIS - tests compositional difference between
groups “The ANOSIM statistic “R” compares the mean of ranked
dissimilarities between groups to the mean of ranked dissimilarities
within groups. An R value close to “1.0” suggests dissimilarity between
groups while an R value close to “0” suggests an even distribution of
high and low ranks within and between groups” (GUSTAME). In other words,
the higher the R value, the more dissimilar your groups are in terms of
microbial community composition.

### Kruskal-Wallis testing for significant Taxa

### Spearman corrs of taxa with CH4pw

## Taxa faceted PCoA Plot with Bray-curtis distances

### Methano-analysis

### Subset phyloseq object for Methanogens

## Taxa bar plots (June - Methanogenic Orders)

### Set up data to visualise relative abundance and group phyla

### Plot

#### Taxa Stacked METHANO BY GENUS

### Alpha Diversity

#### now for mgens

\#Barplot Counts

\#Prep PS object for SpeicEasi

# Core & shared taxa
