# humann3 analysis TFG

library(ggplot2)
library(vegan)
require(gplots)
library(pheatmap)
library(dplyr)
library(viridis)

setwd("/media/sequentia/visitors/visitor8/TFG")

# read files

PD_genefamilies <- read.table("read-based/Functional_classification/humann3/PD_genefamilies.tsv", 
                              stringsAsFactors = F, sep="\t", quote="", header = TRUE)

PD_pathabundance <- read.table("read-based/Functional_classification/humann3/PD_pathabundance.tsv", 
                               stringsAsFactors = F, sep="\t", quote="", header = TRUE)

control_genefamilies <- read.table("read-based/Functional_classification/humann3/CONTROL_genefamilies.tsv", 
                                   stringsAsFactors = F, sep="\t", quote="", header = TRUE)

control_pathabundance <- read.table("read-based/Functional_classification/humann3/CONTROL_pathabundance.tsv", 
                               stringsAsFactors = F, sep="\t", quote="", header = TRUE)


## Qualitative study

#########################
# heatmap pathabundance #
#########################

PD_path <- PD_pathabundance
control_path <- control_pathabundance

colnames(PD_path) <- c("Pathway","SRR15853012", "SRR15852941", "SRR15852948","SRR15853009", "SRR15853002",
                       "SRR15852986", "SRR15852987", "SRR15852998", "SRR15853011","SRR15853000",
                       "SRR15853017","SRR15853005", "SRR15852988", "SRR15852994","SRR15852995", 
                       "SRR15852997", "SRR15852989","SRR15852990","SRR15853001","SRR15852992")

colnames(control_path) <- c("Pathway", "SRR15852951", "SRR15852922", "SRR15852955", "SRR15852956","SRR15852962",
                            "SRR15852963", "SRR15852982","SRR15853026","SRR15853027","SRR15852932",
                            "SRR15852950", "SRR15852952","SRR15852954","SRR15852957","SRR15852960",
                            "SRR15852993", "SRR15852934","SRR15852968","SRR15852969","SRR15853028")

# remove rows having UNINTEGRATED
PD_removed <- PD_path[!grepl("UNINTEGRATED", PD_path$Pathway),]
control_removed <- control_path[!grepl("UNINTEGRATED", control_path$Pathway),]

# check for host contamination (Homo sapiens) in the sample (no contamination found)
mm<-apply(control_removed,1,paste,collapse=" ")
grep("s__homo", mm)
control_removed[2089,]

# remove vertical line (removing the pathways classified for each specie)
PD_removed <- PD_removed[grep("\\|",PD_removed$Pathway, invert = T),] # invert: remove rows having that pattern
control_removed <- control_removed[grep("\\|",control_removed$Pathway, invert = T),]

# merge both files
merged <- merge(PD_removed, control_removed, by = "Pathway", all = T )
merged <- merge(PD_removed, control_removed, by = "Pathway", all = T )
merged2 <- merge(PD_removed, control_removed, by = "Pathway", all = T )


# convert path names into row names, so the data matrix contains only sequence count data
row.names(PD_removed) <- PD_removed$Pathway
PD_path_data <- PD_removed
PD_path_data <- PD_path_data[,-1]

row.names(control_removed) <- control_removed$Pathway
control_path_data <- control_removed
control_path_data <- control_path_data[,-1]

row.names(merged) <- merged$Pathway
merged <- merged[,-1]
# replace NA by 0's
merged[is.na(merged)] <- 0

functional_path_matrix <- merged

# plot
pheatmap(PD_path_data, scale = 'row',show_rownames = F, main = "PD pathway abundance")
pheatmap(control_path_data, scale = 'row',show_rownames = F, main = "control pathway abundance")
pheatmap(merged, scale = 'row',show_rownames = F, main = "Pathway abundance")



###########################
# z-score standardization #
###########################

data_normalized <- merged2
# normalize only numeric values
data_normalized <- scale(data_normalized[,2:41])

# extract the pathway name column
column <- merged2[,1, drop = FALSE]

# bind pathway column and normalized values data.frame
data_norm <- cbind(column, data_normalized)

# replace NA by 0's
data_norm[is.na(data_norm)] <- 0

# reshape data into long format
data_long_norm <-  melt(data_norm, id.vars = "Pathway", variable.name = "sample")

## Add samples status

# convert sample column factor into characters
dataFac <- data_long_norm
dataFac <- data.frame(lapply(dataFac, as.character), stringsAsFactors = FALSE)
# adding column status based on other column
dataFac <- dataFac %>%
  mutate(status = case_when(
    endsWith(sample, "3012") ~ "PD",
    endsWith(sample, "2941") ~ "PD",
    endsWith(sample, "2948") ~ "PD",
    endsWith(sample, "3009") ~ "PD",
    endsWith(sample, "3002") ~ "PD",
    endsWith(sample, "2986") ~ "PD",
    endsWith(sample, "2987") ~ "PD",
    endsWith(sample, "2998") ~ "PD",
    endsWith(sample, "3011") ~ "PD",
    endsWith(sample, "3000") ~ "PD",
    endsWith(sample, "3017") ~ "PD",
    endsWith(sample, "3005") ~ "PD",
    endsWith(sample, "2988") ~ "PD",
    endsWith(sample, "2994") ~ "PD",
    endsWith(sample, "2995") ~ "PD",
    endsWith(sample, "2997") ~ "PD",
    endsWith(sample, "2989") ~ "PD",
    endsWith(sample, "2990") ~ "PD",
    endsWith(sample, "3001") ~ "PD",
    endsWith(sample, "2992") ~ "PD",
    endsWith(sample, "2951") ~ "control",
    endsWith(sample, "2922") ~ "control",
    endsWith(sample, "2955") ~ "control",
    endsWith(sample, "2956") ~ "control",
    endsWith(sample, "2962") ~ "control",
    endsWith(sample, "2963") ~ "control",
    endsWith(sample, "2982") ~ "control",
    endsWith(sample, "3026") ~ "control",
    endsWith(sample, "3027") ~ "control",
    endsWith(sample, "2932") ~ "control",
    endsWith(sample, "2950") ~ "control",
    endsWith(sample, "2952") ~ "control",
    endsWith(sample, "2954") ~ "control",
    endsWith(sample, "2957") ~ "control",
    endsWith(sample, "2960") ~ "control",
    endsWith(sample, "2993") ~ "control",
    endsWith(sample, "2934") ~ "control",
    endsWith(sample, "2968") ~ "control",
    endsWith(sample, "2969") ~ "control",
    endsWith(sample, "3028") ~ "control"
  ))
View(dataFac)

# from character to factor

dataFin <- dataFac
dataFin$Pathway <- as.factor(dataFin$Pathway)
dataFin$sample <- as.factor(dataFin$sample)

# from character to numeric
dataFin$value <- as.numeric(dataFin$value)

# plot
p <- ggplot(dataFin, aes(sample,Pathway, fill =value)) +geom_tile() + xlab(label = "Samples") + ylab("") + 
  theme(axis.text.x = element_text(angle =90 )) +facet_grid(.~ status, scales = "free_x", space = "free_x") + scale_fill_viridis()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
p

# plot 2

# make a data frame with 2 columns, first is the samples and second is the status
my_sample_col2 <- data.frame(status = rep(c("PD", "control"), c(20,20)))
my_sample_col <- data.frame(samples = rep(c("SRR15853012", "SRR15852941", "SRR15852948","SRR15853009", "SRR15853002",
                                            "SRR15852986", "SRR15852987", "SRR15852998", "SRR15853011","SRR15853000",
                                            "SRR15853017","SRR15853005", "SRR15852988", "SRR15852994","SRR15852995", 
                                            "SRR15852997", "SRR15852989","SRR15852990","SRR15853001","SRR15852992",
                                            "SRR15852951", "SRR15852922", "SRR15852955", "SRR15852956","SRR15852962",
                                            "SRR15852963", "SRR15852982","SRR15853026","SRR15853027","SRR15852932",
                                            "SRR15852950", "SRR15852952","SRR15852954","SRR15852957","SRR15852960",
                                            "SRR15852993", "SRR15852934","SRR15852968","SRR15852969","SRR15853028")))
my_col <- cbind(my_sample_col,my_sample_col2)

row.names(my_col) <- my_col$samples

my_col$samples <- NULL

pheatmap(merged, annotation_col = my_col, scale = 'row',show_rownames = F, width=15, height=8,main = "Pathway abundance")


tiff("test.tiff2", units="in", width=15, height=8, res=300)
pheatmap(merged, annotation_col = my_col, scale = 'row',show_rownames = F, main = "Pathway abundance")
dev.off()



#######################
# heapmap GENE_FAMILY #
#######################

PD_gene <- PD_genefamilies
control_gene <- control_genefamilies

colnames(PD_gene) <- c("Gene_family", "SRR15853012", "SRR15852941", "SRR15852948","SRR15853009", "SRR15853002",
                       "SRR15852986", "SRR15852987", "SRR15852998", "SRR15853011","SRR15853000",
                       "SRR15853017","SRR15853005", "SRR15852988", "SRR15852994","SRR15852995", 
                       "SRR15852997", "SRR15852989","SRR15852990","SRR15853001","SRR15852992")

colnames(control_gene) <- c("Gene_family", "SRR15852951", "SRR15852922", "SRR15852955", "SRR15852956","SRR15852962",
                            "SRR15852963", "SRR15852982","SRR15853026","SRR15853027","SRR15852932",
                            "SRR15852950", "SRR15852952","SRR15852954","SRR15852957","SRR15852960",
                            "SRR15852993", "SRR15852934","SRR15852968","SRR15852969","SRR15853028")


# remove vertical line (meaning removing the pathways classified for each specie)
PD_gene <- PD_gene[grep("\\|",PD_gene$Gene_family, invert=T),] # invert: remove rows having that pattern
control_gene <- control_gene[grep("\\|",control_gene$Gene_family, invert=T),] # invert: remove rows having that pattern



## merge both files
merged_genes <- merge(PD_gene, control_gene, by = "Gene_family", all = T)

# convert gene names into row names, so the data matrix contains only sequence count data
row.names(merged_genes) <- merged_genes$Gene_family
merged_genes <- merged_genes[,2:41]
# replace NA by 0's
merged_genes[is.na(merged_genes)] <- 0



#### (go down to the R file)

# Genes as row names
rownames(PD_gene) <- PD_gene$Gene_family 
PD_gene <- PD_gene[,2:21]

rownames(control_gene) <- control_gene$Gene_family 
control_gene <- control_gene[,2:21]


# compute abundances
PD_gene_percent <- t(t(PD_gene)/colSums(PD_gene))*100
control_gene_percent <- t(t(control_gene)/colSums(control_gene))*100
general_gene_percent <- t(t(merged_genes)/colSums(merged_genes))*100
# max(PD_gene_percent) : 1.681524
# max(control_gene_percent) : 0.8173253
# max(general_gene_percent) : 1.681524

# remove species with low abundance, below 0.005 and 0.05 for the merged one
PD_gene_perc_filt <- PD_gene_percent[!apply(PD_gene_percent, 1, function(x){sum(x<0.005)}) == 20,] # 20 is num of samples
control_gene_perc_filt <- control_gene_percent[!apply(control_gene_percent, 1, function(x){sum(x<0.005)}) == 20,] # 20 is num of samples
general_gene_perc_filt <- general_gene_percent[!apply(general_gene_percent, 1, function(x){sum(x<0.05)}) == 40,] # 40 is num of samples

# summarize % of removed gene family
PD_gene_perc_filt <- rbind(PD_gene_perc_filt, 100-colSums(PD_gene_perc_filt))
rownames(PD_gene_perc_filt)[nrow(PD_gene_perc_filt)] <- "Other"

control_gene_perc_filt <- rbind(control_gene_perc_filt, 100-colSums(control_gene_perc_filt))
rownames(control_gene_perc_filt)[nrow(control_gene_perc_filt)] <- "Other"

general_gene_perc_filt <- rbind(general_gene_perc_filt, 100-colSums(general_gene_perc_filt))
rownames(general_gene_perc_filt)[nrow(general_gene_perc_filt)] <- "Other"


# plot
pheatmap(PD_gene_perc_filt, scale = 'row',show_rownames = F, main = " PD Gene family abundance")
pheatmap(control_gene_perc_filt, scale = 'row',show_rownames = F, main = " control Gene family abundance")


# plot 2

# make a data frame with 2 columns, first is the samples and second is the status
my_sample_col2 <- data.frame(status = rep(c("PD", "control"), c(20,20)))
my_sample_col <- data.frame(samples = rep(c("SRR15853012", "SRR15852941", "SRR15852948","SRR15853009", "SRR15853002",
                                            "SRR15852986", "SRR15852987", "SRR15852998", "SRR15853011","SRR15853000",
                                            "SRR15853017","SRR15853005", "SRR15852988", "SRR15852994","SRR15852995", 
                                            "SRR15852997", "SRR15852989","SRR15852990","SRR15853001","SRR15852992",
                                            "SRR15852951", "SRR15852922", "SRR15852955", "SRR15852956","SRR15852962",
                                            "SRR15852963", "SRR15852982","SRR15853026","SRR15853027","SRR15852932",
                                            "SRR15852950", "SRR15852952","SRR15852954","SRR15852957","SRR15852960",
                                            "SRR15852993", "SRR15852934","SRR15852968","SRR15852969","SRR15853028")))
my_col <- cbind(my_sample_col,my_sample_col2)

row.names(my_col) <- my_col$samples

my_col$samples <- NULL

pheatmap(general_gene_perc_filt, annotation_col = my_col, scale = 'row',show_rownames = F,fontsize = 12, main = "Genes family abundance", res=300)

tiff("test.tiff", units="in", width=15, height=8, res=300)
pheatmap(general_gene_perc_filt, annotation_col = my_col, scale = 'row',show_rownames = F,fontsize = 12, main = "Genes family abundance", res=300)
dev.off()


#####################
# ML matrix
#####################

# normalize the data into a scale of 0 to 100 (computing the abundance)
functional_path_matrix <- t(t(functional_path_matrix)/colSums(functional_path_matrix))*100

# check that columns sum up to 100
colSums(functional_path_matrix[,1:40])

## metadata

# To aggregate the metadata and to be allow to use the machine learning algorithmns we need
# to have our samples as rows and the species name as columns
# so we have to transpose our matrix

# Transpose matrix
funct_matrix_t <- data.frame(t(functional_path_matrix))

# metadata vector
condition <- c("PD","PD","PD","PD","PD","PD","PD","PD","PD","PD","PD","PD","PD","PD","PD",
                   "PD","PD","PD","PD","PD","control","control","control","control","control",
                   "control","control","control","control","control","control","control",
                   "control","control","control","control","control","control","control","control")

condition2 <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

# Bind metadata and data.frame
funct_matrix_t <- cbind(funct_matrix_t,condition2)

# take a look at the column condition
grep("condition2", colnames(funct_matrix_t))
funct_matrix_t[,486:487]
class(funct_matrix_t)

# Download data frame to CSV file
write.csv(funct_matrix_t, file = '/media/sequentia/visitors/visitor8/TFG/pathway_matrix2.csv')






## Gene family matrix ML (try again to download)

geneFamily_matrix <- merged_genes

# normalize the data into a scale of 0 to 100 (computing the abundance)
geneFamily_matrix <- t(t(geneFamily_matrix)/colSums(geneFamily_matrix))*100

# check that columns sum up to 100
colSums(geneFamily_matrix[,1:40])

## metadata

# To aggregate the metadata and to be allow to use the machine learning algorithms we need
# to have our samples as rows and the species name as columns
# so we have to transpose our matrix

# Transpose matrix
geneFamily_matrix_t <- data.frame(t(geneFamily_matrix))

# metadata vector
condition <- c("PD","PD","PD","PD","PD","PD","PD","PD","PD","PD","PD","PD","PD","PD","PD",
               "PD","PD","PD","PD","PD","control","control","control","control","control",
               "control","control","control","control","control","control","control",
               "control","control","control","control","control","control","control","control")

# Bind metadata and data.frame
geneFamily_matrix_t <- cbind(geneFamily_matrix_t,condition)

# take a look at the column condition
grep("condition", colnames(geneFamily_matrix_t))
geneFamily_matrix_t[,486:487]
class(geneFamily_matrix_t)

# Download data frame to CSV file
write.csv(geneFamily_matrix_t, file = '/media/sequentia/visitors/visitor8/TFG/geneFamily2_matrix.csv')
