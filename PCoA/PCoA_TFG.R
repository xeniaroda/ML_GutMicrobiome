
# PCoA 


library(phyloseq)
library(ggplot2)
library(vegan)
library(ape)
library(plyr)
library(microbiome)
library(ggpubr)
library("DESeq2")
library(ggprism)

######################
## Taxonomic matrix ##
######################
setwd("/media/sequentia/visitors/visitor8/TFG")


# we load the merged file data, done in centrifuge_TFG.r
# this are the count for each species
load("merged.Rda")


## Diversity analysis
# create phyloseq object
species_phyloseq <- otu_table(merged, taxa_are_rows = T)

# rarefraction curve
rarecurve(t(species_phyloseq), step = 10000, label = T, ylab= "Species")

#####################
# diversity metrics #
#####################

# make some plots to visualize alpha-diversity
rich <- plot_richness(species_phyloseq)

# get alpha-diversity numbers
alpha <- estimate_richness(species_phyloseq)

# check Shannon diversity
alpha$conditions <- c("PD", "PD", "PD", "PD", "PD", "PD", "PD", "PD", "PD", "PD", "PD", "PD"
                      , "PD", "PD", "PD", "PD", "PD", "PD", "PD", "PD","control","control"
                      ,"control","control","control","control","control","control","control"
                      ,"control","control","control","control","control","control","control"
                      ,"control","control","control","control")


ggplot(alpha, aes(x= conditions, y = Shannon))+
  geom_boxplot()+ stat_compare_means(method = "wilcox.test")



###################################
# differential abundance analysis #
###################################

# create phyloseq metadata object
conditions_pyseq <- sample_data(data.frame(row.names=c("SRR15853012", "SRR15852941", "SRR15852948","SRR15853009", "SRR15853002",
                                                      "SRR15852986", "SRR15852987", "SRR15852998", "SRR15853011","SRR15853000",
                                                      "SRR15853017","SRR15853005", "SRR15852988", "SRR15852994","SRR15852995", 
                                                      "SRR15852997", "SRR15852989","SRR15852990","SRR15853001","SRR15852992",
                                                      "SRR15852951", "SRR15852922", "SRR15852955", "SRR15852956","SRR15852962",
                                                      "SRR15852963", "SRR15852982","SRR15853026","SRR15853027","SRR15852932",
                                                      "SRR15852950", "SRR15852952","SRR15852954","SRR15852957","SRR15852960",
                                                      "SRR15852993", "SRR15852934","SRR15852968","SRR15852969","SRR15853028"),
                                           condition=c("PD", "PD", "PD", "PD", "PD", "PD", "PD", "PD", "PD", "PD", "PD", "PD"
                                                       , "PD", "PD", "PD", "PD", "PD", "PD", "PD", "PD","control","control"
                                                       ,"control","control","control","control","control","control","control"
                                                       ,"control","control","control","control","control","control","control"
                                                       ,"control","control","control","control")))
                                           

# make condition as reference
conditions_pyseq$condition <- factor(conditions_pyseq$condition, levels=c("control", "PD"))

pyseq <- phyloseq(conditions_pyseq, species_phyloseq)
dds <- phyloseq_to_deseq2(pyseq, ~ condition)
dds <- DESeq(dds) # consistent with our analysis

# results of the differential abundance analysis
res <- as.data.frame(results(dds))
sig <- res[which(res$padj < 0.05), ]
sig <- as.data.frame(sig)

# permanova to assess significance of PCoA
bray <- phyloseq::distance(pyseq, "bray")
bray <- as.matrix(bray)

p.value.cal <-adonis2(bray ~ condition, data = data.frame(pyseq@sam_data))
p.value <- p.value$`Pr(>F)`[1] # p_value: 0.53

########
# PCoA #
########

pyseq <- phyloseq(conditions_pyseq, otu_table(t(t(species_phyloseq)/colSums(species_phyloseq))*100), TRUE)

# we add a "dummy variable" in sample data
sample_data(pyseq)[,2] <- sample_data(pyseq)[,1]

# plot 1: bray-Curtis beta diversity metric
data.ord <- ordinate(pyseq, "PCoA", "bray")
p1<- plot_ordination(pyseq, data.ord, type = "samples", color = "condition", shape= "condition",
                     title = "Principal Coordinates (PCoA)")+ annotate("text", x=0.32, y=-0.60, label= "p_value = 0.53") + stat_ellipse() 

p1.2 <- plot_ordination(pyseq, data.ord, type = "samples", color = "condition", shape= "condition",
                        title = "Principal Coordinates (PCoA)")+ facet_wrap(~condition,3)


# plot 2: euclidean metric
data.ord.2 <- ordinate(pyseq, "PCoA", "euclidean")
p2<- plot_ordination(pyseq, data.ord.2, type = "samples", color = "condition", shape= "condition",
                     title = "Principal Coordinates (PCoA)")

p2 + facet_wrap(~condition,3)

# plot 3: jaccard metric
data.ord.3 <- ordinate(pyseq, "PCoA", "jaccard")
p3<- plot_ordination(pyseq, data.ord.3, type = "samples", color = "condition", shape= "condition",
                     title = "Principal Coordinates (PCoA)")


p3 + facet_wrap(~condition,3)


# plot 4 : manhattan metric
data.ord.4 <- ordinate(pyseq, "PCoA", "manhattan")
p4<- plot_ordination(pyseq, data.ord.4, type = "samples", color = "condition", shape= "condition",
                     title = "Principal Coordinates (PCoA)")


p4 + facet_wrap(~condition,3)



# CONCLUSION:
# After performing our PCoA with different beta diversity metrics,
# we see that samples from both conditions are overlapping together,
# but still they have their own trend to aggregate separately.
# This means that with basic statistical methods we can not visualize any pattern in our data
# thus, we have more reasons to use AI techniques, which in this way we can discover hidden patterns.


### method 2 (without phyloseq)
# without transforming the data into phyloseq object
dist.matrix <- vegdist(taxo_matrix.num, method = "bray")


# PCoA is not included in vegan, use ape package instead
PCOA <- pcoa(dist.matrix, correction = "cailliez")

# plot the results
biplot.pcoa(PCOA)

biplot.pcoa(PCOA, taxo_matrix.num)


