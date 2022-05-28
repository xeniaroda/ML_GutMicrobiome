
# LOOCV

library(Boruta)
library(caret)
library(neuralnet)
library(devtools)
library(nnet)
library(NeuralNetTools)
library(data.table)
library(randomForest)
library(rpart)
library(rpart.plot)
library(naivebayes)
library(dplyr)
library(psych)
library(e1071)
library(caTools)
library(class)
library(kernlab)
library(C50)
library(klaR)
library(pROC)
library(MLeval)
library(tidyverse)
library(Boruta)
library(caret)
library(randomForest)

setwd("/media/sequentia/visitors/visitor8/TFG")

# possible methods : 
names(getModelInfo())

###########################
# 1. Taxonomic matrix #####
###########################

# read the matrix
taxoCV_matrix <- read.csv("taxonomic_matrix2.csv")

######################
## Data preparation ##
######################

# convert sample names in row names
rownames(taxoCV_matrix) <- taxoCV_matrix$X
taxoCV_matrix$X <- NULL

# factor
taxoCV_matrix$condition2<- factor(taxoCV_matrix$condition2, levels = c(0,1))
levels(taxoCV_matrix$condition2) <- c("control", "PD")


##########
## k-NN ##
##########

set.seed(123)

trcontrol_tax_knn1 <- trainControl(method = "LOOCV",number = 5,savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)

fit_tax_knn1 <- train(condition2 ~., method = "knn", trControl = trcontrol_tax_knn1, tuneLength = 0,
                      metric = "Accuracy", data = taxoCV_matrix)
var <- varImp(fit_tax_knn1, scale = F)
plot <- ggplot2::ggplot(var, top = 20) + xlab("Species")



# AUC

knn.taxo.roc <- evalm(fit_tax_knn1) # 0.4
knn.taxo.roc$roc
knn.taxo.roc$stdres

# unfold (by cv)

v1 <- trainControl(method = "cv", number = 5,savePredictions = TRUE, verboseIter = TRUE)

v  <- train(condition2 ~., method = "knn", trControl = v1, tuneLength = 0, data = taxoCV_matrix)


pred <- v$pred
pred$equal <- ifelse(pred$pred == pred$obs, 1,0)
eachfold1 <- pred %>%
  group_by(Resample) %>%
  summarise_at(vars(equal),
               list(Accuracy = mean))


ggplot(data=eachfold1, aes(x=factor(Resample), y=Accuracy, group = 1)) +
  geom_boxplot(color="maroon") +
  geom_point() +
  theme_minimal()



#################
## Naive Bayes ## 
#################


grep("condition2", colnames(taxoCV_matrix)) # column 5658

trcontrol_tax_nb <- trainControl(method = "LOOCV",savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
fit_tax_nb <- train(condition2 ~., method = "naive_bayes", trControl = trcontrol_tax_nb,
                      metric = "Accuracy", data = taxoCV_matrix)
var2 <-varImp(fit_tax_nb, scale = F)
ggplot2::ggplot(var2, top = 20) + xlab("Species")



# AUC

nb.taxo.roc <- evalm(fit_tax_nb)
nb.taxo.roc$roc
nb.taxo.roc$stdres

# unfold

c1 <- trainControl(method = "cv",savePredictions = TRUE, number = 5,verboseIter = TRUE)
c <- train(condition2 ~., method = "naive_bayes", trControl = c1,
                    metric = "Accuracy",tuneLength = 0, data = taxoCV_matrix)

pred <- c$pred
pred$equal <- ifelse(pred$pred == pred$obs, 1,0)
eachfold2 <- pred %>%
  group_by(Resample) %>%
  summarise_at(vars(equal),
               list(Accuracy = mean))

ggplot(data=eachfold2, aes(x=factor(Resample), y=Accuracy, group = 1)) +
  geom_boxplot(color="maroon") +
  geom_point() +
  theme_minimal()



################################
## Support vector machine SVM ##
################################

trcontrol_tax_SVM <- trainControl(method = "LOOCV",savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
fit_tax_SVM <- train(condition2 ~., method = "svmLinear2", trControl = trcontrol_tax_SVM,
                    metric = "Accuracy", data = taxoCV_matrix)

var3 <-varImp(fit_tax_SVM, scale = F)
ggplot2::ggplot(var3, top = 20) + xlab("Species")

# AUC
SVM.taxo.roc <- evalm(fit_tax_SVM)
SVM.taxo.roc$roc
SVM.taxo.roc$stdres

# unfold

t1 <- trainControl(method = "cv",savePredictions = TRUE, number = 5 ,classProbs = TRUE, verboseIter = TRUE)
t <- train(condition2 ~., method = "svmLinear2", trControl = t1, tuneLength = 0,
                     metric = "Accuracy", data = taxoCV_matrix)


pred <- t$pred
pred$equal <- ifelse(pred$pred == pred$obs, 1,0)
eachfold3 <- pred %>%
  group_by(Resample) %>%
  summarise_at(vars(equal),
               list(Accuracy = mean))

ggplot(data=eachfold3, aes(x=factor(Resample), y=Accuracy, group = 1)) +
  geom_boxplot(color="maroon") +
  geom_point() +
  theme_minimal()

###################
## Decision tree ##
###################

set.seed(123)
trcontrol_tax_dectree <- trainControl(method = "LOOCV", number = 5, savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
fit_tax_dectree <- train(condition2 ~., method = "C5.0", trControl = trcontrol_tax_dectree,
                     metric = "usage", data = taxoCV_matrix)

var4 <-varImp(fit_tax_dectree, scale = F)
ggplot2::ggplot(var4, top = 20) + xlab("Species")

# AUC
dec.taxo.roc <- evalm(fit_tax_dectree)
dec.taxo.roc$roc
dec.taxo.roc$stdres


# unfold
q1<- trainControl(method = "cv", number = 5, savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
q <- train(condition2 ~., method = "C5.0Tree", trControl = q1, tuneLength = 0,
                         data = taxoCV_matrix)

pred <- q$pred
pred$equal <- ifelse(pred$pred == pred$obs, 1,0)
eachfold4 <- pred %>%
  group_by(Resample) %>%
  summarise_at(vars(equal),
               list(Accuracy = mean))

ggplot(data=eachfold4, aes(x=factor(Resample), y=Accuracy, group = 1)) +
  geom_boxplot(color="maroon") +
  geom_point() +
  theme_minimal()

####################
## Random forest ##
###################
set.seed(123)
trcontrol_tax_rf <- trainControl(method = "LOOCV", savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
fit_tax_rf2 <- train(condition2 ~., method = "rf", trControl = trcontrol_tax_rf,
                    metric = "Accuracy", data = taxoCV_matrix, imp = 1)

var5 <- varImp(fit_tax_rf2, scale = F)
ggplot2::ggplot(var5, top = 20) + xlab("Species")

# AUC
rf.taxo.roc <- evalm(fit_tax_rf2)
rf.taxo.roc$roc
rf.taxo.roc$stdres

# unfold
x <- taxoCV_matrix[,1:5657]
mtry <- sqrt(ncol(x))
tunegrid <- expand.grid(.mtry=mtry)
r1 <- trainControl(method = "cv", number = 5, savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
r <- train(condition2 ~., method = "rf", trControl = r1, tuneLength = 0,
                     metric = "Accuracy", data = taxoCV_matrix, tuneGrid = tunegrid)

pred <- r$pred
pred$equal <- ifelse(pred$pred == pred$obs, 1,0)
eachfold5 <- pred %>%
  group_by(Resample) %>%
  summarise_at(vars(equal),
               list(Accuracy = mean))

ggplot(data=eachfold5, aes(x=factor(Resample), y=Accuracy, group = 1)) +
  geom_boxplot(color="maroon") +
  geom_point() +
  theme_minimal()


########################
# 2. Pathway matrix ####
########################

# possible methods : 
names(getModelInfo())

# read the matrix
pathCV_matrix <- read.csv("pathway_matrix2.csv")

# factor
pathCV_matrix$condition2<- factor(pathCV_matrix$condition2, levels = c(0,1))
levels(pathCV_matrix$condition2) <- c("control", "PD")

######################
## Data preparation ##
######################

# convert sample names in row names
rownames(pathCV_matrix) <- pathCV_matrix$X
pathCV_matrix$X <- NULL



##########
## k-NN ##
##########

set.seed(123)

trcontrol_path_knn1 <- trainControl(method = "LOOCV", savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
fit_path_knn1 <- train(condition2 ~., method = "knn", trControl = trcontrol_path_knn1,
                      metric = "Accuracy", data = pathCV_matrix)


var6 <-varImp(fit_path_knn1, scale = F)

# AUC
knn.path.roc <- evalm(fit_path_knn1)
knn.path.roc$stdres

#################
## Naive Bayes ##
#################


grep("condition2", colnames(pathCV_matrix)) # column 487

trcontrol_path_nb <- trainControl(method = "LOOCV", savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
fit_path_nb <- train(condition2 ~., method = "naive_bayes", trControl = trcontrol_path_nb,
                    metric = "Accuracy", data = pathCV_matrix)

var7 <-varImp(fit_path_nb, scale = F)
# AUC
nb.path.roc <- evalm(fit_path_nb)
nb.taxo.roc$stdres


################################
## Support vector machine SVM ##
################################

set.seed(123)
trcontrol_path_SVM <- trainControl(method = "LOOCV", savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
fit_path_SVM <- train(condition2 ~., method = "svmLinear2", trControl = trcontrol_path_SVM,
                     metric = "Accuracy", data = pathCV_matrix)

var8 <-varImp(fit_path_SVM, scale = F)
# AUC
SVM.path.roc <- evalm(fit_path_SVM)
SVM.path.roc$stdres

###################
## Decision tree ##
###################

set.seed(123)
trcontrol_path_dectree <- trainControl(method = "LOOCV", savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
fit_path_dectree <- train(condition2 ~., method = "C5.0", trControl = trcontrol_path_dectree,
                         metric = "Accuracy", data = pathCV_matrix)

var9 <-varImp(fit_path_dectree, scale = F)
# AUC
dec.path.roc <- evalm(fit_path_dectree)
dec.path.roc$stdres


####################
## Random forest ##
###################
set.seed(123)
trcontrol_path_rf <- trainControl(method = "LOOCV", savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
fit_path_rf2 <- train(condition2 ~., method = "rf", trControl = trcontrol_path_rf,
                     metric = "Accuracy", data = pathCV_matrix)

var10 <-varImp(fit_path_rf2, scale = F)

# AUC
rf.path.roc <- evalm(fit_path_rf2)
rf.path.roc$stdres


#####################
# 3. Gene family ####
#####################

# possible methods : 
names(getModelInfo())

# read the matrix
geneCV_matrix <- read.csv("geneFamilySmall2_matrix.csv")

# factor
geneCV_matrix$condition2<- factor(geneCV_matrix$condition2, levels = c(0,1))
levels(geneCV_matrix$condition2) <- c("control", "PD")

# convert sample names in row names
rownames(geneCV_matrix) <- geneCV_matrix$X
geneCV_matrix$X <- NULL



##########
## k-NN ##
##########

set.seed(123)

trcontrol_gene_knn1 <- trainControl(method = "LOOCV",  savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
fit_gene_knn1 <- train(condition2 ~., method = "knn", trControl = trcontrol_gene_knn1,
                       metric = "Accuracy", data = geneCV_matrix)

var11 <- varImp(fit_gene_knn1, scale = F)

# AUC
knn.gene.roc <- evalm(fit_gene_knn1)
knn.gene.roc$stdres

#################
## Naive Bayes ##
#################

trcontrol_gene_nb <- trainControl(method = "LOOCV", savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
fit_gene_nb <- train(condition2 ~., method = "naive_bayes", trControl = trcontrol_gene_nb,
                     metric = "Accuracy", data = geneCV_matrix)

var12 <- varImp(fit_gene_knn1, scale = F)

# AUC
nb.gene.roc <- evalm(fit_gene_nb)
nb.gene.roc$stdres


################################
## Support vector machine SVM ##
################################
set.seed(123)
trcontrol_gene_SVM <- trainControl(method = "LOOCV",savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
fit_gene_SVM <- train(condition2 ~., method = "svmLinear2", trControl = trcontrol_gene_SVM,
                      metric = "Accuracy", data = geneCV_matrix)

var13 <- varImp(fit_gene_SVM, scale = F)

# AUC
SVM.gene.roc <- evalm(fit_gene_SVM)
SVM.gene.roc$stdres

###################
## Decision tree ##
###################

set.seed(123)
trcontrol_gene_dectree <- trainControl(method = "LOOCV", savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
fit_gene_dectree <- train(condition2 ~., method = "C5.0", trControl = trcontrol_gene_dectree,
                          metric = "Accuracy", data = geneCV_matrix)

var14 <- varImp(fit_gene_dectree, scale = F)

# AUC
dec.gene.roc <- evalm(fit_gene_dectree)
dec.gene.roc$stdres


####################
## Random forest ##
###################
set.seed(123)
trcontrol_gene_rf <- trainControl(method = "LOOCV", savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
fit_gene_rf2 <- train(condition2 ~., method = "rf", trControl = trcontrol_gene_rf,
                      metric = "Accuracy", data = geneCV_matrix)

var15 <- varImp(fit_gene_rf2, scale = F)

# AUC
rf.gene.roc <- evalm(fit_gene_rf2)
rf.gene.roc$stdres


################ SECOND PART ####################

########################
## Feature selection ##
#######################
# Boruta algorithm

### taxonomic
set.seed(111)
boruta.taxo <- Boruta(condition2 ~., data = taxoCV_matrix, doTrace = 2, maxRuns = 500)

# list of confirmed attributes
getSelectedAttributes(boruta.taxo, withTentative = F)

# dataframe of results
boruta.df <- attStats(boruta.taxo)
print(boruta.df)

# Data partition
set.seed(222)
ind <- sample(2, nrow(taxoCV_matrix), replace = T, prob = c(0.6, 0.4))
train.taxo <- taxoCV_matrix[ind==1,] # 26 5658
test.taxo <- taxoCV_matrix[ind==2,] # 14 5658

set.seed(333)
rf.taxo1 <- trainControl(method = "LOOCV", savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
rf.taxo11 <- train(condition2 ~., method = "rf", trControl = trcontrol_tax_rf,
                     metric = "Accuracy", data = train.taxo) #acc = 0.4615385

rf.feature.test <- train(condition2 ~Arthrobacter.phage.KellEzio +
                      Escherichia.phage.PBECO.4 + Flexistipes.sinusarabici + Lactococcus.phage.50101 +
                      Salmonella.phage.vB_SosS_Oslo + Cellulophaga.phage.phi38.1 + Escherichia.phage.vB_EcoM_Alf5
                    + Geobacillus.sp..JS12 + Mycobacterium.phage.Charlie + Yellowstone.lake.mimivirus, method = "rf", trControl = trcontrol_tax_rf,
                    metric = "Accuracy", data = test.taxo)


### pathway abundance

set.seed(111)
boruta.path <- Boruta(condition2 ~., data = pathCV_matrix, doTrace = 2, maxRuns = 600)

# take decision on tentative attributes
final.boruta.path <- TentativeRoughFix(boruta.path)

# list of confirmed attributes
getSelectedAttributes(final.boruta.path, withTentative = F)


### gene family abundance

set.seed(111)
boruta.gene <- Boruta(condition2 ~., data = geneCV_matrix, doTrace = 2, maxRuns = 600)

# list of confirmed attributes
getSelectedAttributes(boruta.gene, withTentative = F)





##############
## boxplots ##
##############

# Taxonomic

taxo.feature <- taxoCV_matrix
taxo.feature <- t(taxo.feature)

# convert row names into the first column
library(tibble)
taxo.feature <- as.data.frame(taxo.feature)
data.taxo <- tibble::rownames_to_column(taxo.feature, "species")

# reshape data into long format
data.taxo <- melt(data.taxo, id.vars = "species", variable.name = "sample")

data.feature <- data.taxo
data.feature <- data.frame(lapply(data.feature, as.character), stringsAsFactors = FALSE)
# adding column status based on other column
data.feature <- data.feature %>%
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
View(data.feature)  

# from character to factor

data.feature$species <- as.factor(data.feature$species)
data.feature$sample <- as.factor(data.feature$sample)

# from character to numeric
data.feature$value <- as.numeric(data.feature$value)

data.feature2 <- data.feature

# get rows having only the important features
rm(data.feature2)

pat <- paste(c("Arthrobacter.phage.KellEzio", "Escherichia.phage.PBECO.4","Flexistipes.sinusarabici",
               "Lactococcus.phage.50101","Salmonella.phage.vB_SosS_Oslo","Cellulophaga.phage.phi38.1",
               "Escherichia.phage.vB_EcoM_Alf5","Geobacillus.sp..JS12", "Mycobacterium.phage.Charlie",
               "Yellowstone.lake.mimivirus" ), collapse= "|")
data.feature2 <- data.feature2[grepl(pat, data.feature2$species),]

# boxplot
p <-ggplot(data = data.feature2, aes(x = species, y = value)) + geom_boxplot(aes(fill =status ))

p + facet_wrap(~ species, scales = "free")


## pathway

path.feature <- pathCV_matrix
path.feature <- t(path.feature)

# convert row names into the first column
library(tibble)
path.feature <- as.data.frame(path.feature)
data.path <- tibble::rownames_to_column(path.feature, "path")

# reshape data into long format
data.path <- melt(data.path, id.vars = "path", variable.name = "sample")

data.feature.path <- data.path
data.feature.path <- data.frame(lapply(data.feature.path, as.character), stringsAsFactors = FALSE)
# adding column status based on other column
data.feature.path <- data.feature.path %>%
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
View(data.feature.path)  

# from character to factor

data.feature.path$path <- as.factor(data.feature.path$path)
data.feature.path$sample <- as.factor(data.feature.path$sample)

# from character to numeric
data.feature.path$value <- as.numeric(data.feature.path$value)

data.feature2.path <- data.feature.path

# get rows having only the important features

pat2 <- paste(c("PWY.7237..myo...chiro..and.scyllo.inositol.degradation", "PWY0.1221..putrescine.degradation.II",
                "RUMP.PWY..formaldehyde.oxidation.I"), collapse= "|")
data.feature2.path <- data.feature2.path[grepl(pat2, data.feature2.path$path),]



# boxplot
rm(u)
u <- ggplot(data = data.feature2.path, aes(x = path, y = value)) + geom_boxplot(aes(fill =status ))
u + facet_wrap(~ path, scales = "free")



## gene

gene.feature <- geneCV_matrix
gene.feature <- t(gene.feature)

# convert row names into the first column
gene.feature <- as.data.frame(gene.feature)
data.gene <- tibble::rownames_to_column(gene.feature, "gene_family")

# reshape data into long format
data.gene <- melt(data.gene, id.vars = "gene_family", variable.name = "sample")

data.feature.gene <- data.gene
data.feature.gene <- data.frame(lapply(data.feature.gene, as.character), stringsAsFactors = FALSE)
# adding column status based on other column
data.feature.gene <- data.feature.gene %>%
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
View(data.feature.gene)  

# from character to factor

data.feature.gene$gene_family <- as.factor(data.feature.gene$gene_family)
data.feature.gene$sample <- as.factor(data.feature.gene$sample)

# from character to numeric
data.feature.gene$value <- as.numeric(data.feature.gene$value)

data.feature2.gene <- data.feature.gene

# get rows having only the important features

pat3 <- paste(c("UniRef90_A0A078ST01", "UniRef90_A0A139L4S8",
                "UniRef90_E2N7L9"), collapse= "|")
data.feature2.gene <- data.feature2.gene[grepl(pat3, data.feature2.gene$gene_family),]



# boxplot

e <- ggplot(data = data.feature2.gene, aes(x = gene_family, y = value)) + geom_boxplot(aes(fill =status ))
e + facet_wrap(~ gene_family, scales = "free")



###############################
## ML with selected features ##
###############################

# taxonomic 

taxo.feat <- taxoCV_matrix

##########
## k-NN ##
##########

set.seed(123)

trcontrol_tax_knn1.feat <- trainControl(method = "LOOCV",savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)

fit_tax_knn1.feat <- train(condition2 ~
                        Arthrobacter.phage.KellEzio +Escherichia.phage.PBECO.4 +Flexistipes.sinusarabici +
               Lactococcus.phage.50101 +Salmonella.phage.vB_SosS_Oslo +Cellulophaga.phage.phi38.1 +
               Escherichia.phage.vB_EcoM_Alf5 + Geobacillus.sp..JS12 + Mycobacterium.phage.Charlie +
               Yellowstone.lake.mimivirus ,method = "knn", trControl = trcontrol_tax_knn1.feat, tuneLength = 0,
                      metric = "Accuracy", data = taxoCV_matrix)
var <- varImp(fit_tax_knn1.feat, scale = F)


# AUC

knn.taxo.roc.feat <- evalm(fit_tax_knn1.feat) # 0.78
knn.taxo.roc.feat$stdres


#################
## Naive Bayes ## 
#################


trcontrol_tax_nb.feat <- trainControl(method = "LOOCV",savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
fit_tax_nb.feat <- train(condition2 ~ Arthrobacter.phage.KellEzio +Escherichia.phage.PBECO.4 +Flexistipes.sinusarabici +
                           Lactococcus.phage.50101 +Salmonella.phage.vB_SosS_Oslo +Cellulophaga.phage.phi38.1 +
                           Escherichia.phage.vB_EcoM_Alf5 + Geobacillus.sp..JS12 + Mycobacterium.phage.Charlie +
                           Yellowstone.lake.mimivirus, method = "naive_bayes", trControl = trcontrol_tax_nb.feat,
                    metric = "Accuracy", data = taxoCV_matrix)
var2 <-varImp(fit_tax_nb.feat, scale = F)




# AUC

nb.taxo.roc.feat <- evalm(fit_tax_nb.feat)
nb.taxo.roc.feat$stdres


################################
## Support vector machine SVM ##
################################

trcontrol_tax_SVM.feat <- trainControl(method = "LOOCV",savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
fit_tax_SVM.feat <- train(condition2 ~ Arthrobacter.phage.KellEzio +Escherichia.phage.PBECO.4 +Flexistipes.sinusarabici +
                            Lactococcus.phage.50101 +Salmonella.phage.vB_SosS_Oslo +Cellulophaga.phage.phi38.1 +
                            Escherichia.phage.vB_EcoM_Alf5 + Geobacillus.sp..JS12 + Mycobacterium.phage.Charlie +
                            Yellowstone.lake.mimivirus, method = "svmLinear2", trControl = trcontrol_tax_SVM.feat,
                     metric = "Accuracy", data = taxoCV_matrix)

varImp(fit_tax_SVM.feat, scale = F)

# AUC
SVM.taxo.roc.feat <- evalm(fit_tax_SVM.feat)
SVM.taxo.roc.feat$stdres


###################
## Decision tree ##
###################

set.seed(123)
trcontrol_tax_dectree.feat <- trainControl(method = "LOOCV", savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
fit_tax_dectree.feat <- train(condition2 ~ Arthrobacter.phage.KellEzio +Escherichia.phage.PBECO.4 +Flexistipes.sinusarabici +
                                Lactococcus.phage.50101 +Salmonella.phage.vB_SosS_Oslo +Cellulophaga.phage.phi38.1 +
                                Escherichia.phage.vB_EcoM_Alf5 + Geobacillus.sp..JS12 + Mycobacterium.phage.Charlie +
                                Yellowstone.lake.mimivirus, method = "C5.0", trControl = trcontrol_tax_dectree.feat,
                         metric = "usage", data = taxoCV_matrix)

varImp(fit_tax_dectree.feat, scale = F)


# AUC
dec.taxo.roc.feat <- evalm(fit_tax_dectree.feat)
dec.taxo.roc.feat$stdres


####################
## Random forest ##
###################
set.seed(123)
trcontrol_tax_rf.feat <- trainControl(method = "LOOCV", savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
fit_tax_rf2.feat <- train(condition2 ~ Arthrobacter.phage.KellEzio +Escherichia.phage.PBECO.4 +Flexistipes.sinusarabici +
                            Lactococcus.phage.50101 +Salmonella.phage.vB_SosS_Oslo +Cellulophaga.phage.phi38.1 +
                            Escherichia.phage.vB_EcoM_Alf5 + Geobacillus.sp..JS12 + Mycobacterium.phage.Charlie +
                            Yellowstone.lake.mimivirus, method = "rf", trControl = trcontrol_tax_rf.feat,
                     metric = "Accuracy", data = taxoCV_matrix, imp = 1)

varImp(fit_tax_rf2.feat, scale = F)

# AUC
rf.taxo.roc.feat <- evalm(fit_tax_rf2.feat)
rf.taxo.roc.feat$stdres

  
  
  
# pathway

  
##########
## k-NN ##
##########

set.seed(123)

trcontrol_path_knn1.feat <- trainControl(method = "LOOCV", savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
fit_path_knn1.feat <- train(condition2 ~ PWY.7237..myo...chiro..and.scyllo.inositol.degradation + PWY0.1221..putrescine.degradation.II +
                +RUMP.PWY..formaldehyde.oxidation.I, method = "knn", trControl = trcontrol_path_knn1.feat,
                       metric = "Accuracy", data = pathCV_matrix)


varImp(fit_path_knn1.feat, scale = F)

# AUC
knn.path.roc.feat <- evalm(fit_path_knn1.feat)
knn.path.roc.feat$stdres


#################
## Naive Bayes ##
#################

trcontrol_path_nb.feat <- trainControl(method = "LOOCV", savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
fit_path_nb.feat <- train(condition2 ~ PWY.7237..myo...chiro..and.scyllo.inositol.degradation + PWY0.1221..putrescine.degradation.II +
                            +RUMP.PWY..formaldehyde.oxidation.I, method = "naive_bayes", trControl = trcontrol_path_nb.feat,
                     metric = "Accuracy", data = pathCV_matrix)

varImp(fit_path_nb.feat, scale = F)
# AUC
nb.path.roc.feat <- evalm(fit_path_nb.feat)
nb.taxo.roc.feat$stdres


################################
## Support vector machine SVM ##
################################
set.seed(123)
trcontrol_path_SVM.feat <- trainControl(method = "LOOCV", savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
fit_path_SVM.feat <- train(condition2 ~ PWY.7237..myo...chiro..and.scyllo.inositol.degradation + PWY0.1221..putrescine.degradation.II +
                             +RUMP.PWY..formaldehyde.oxidation.I, method = "svmLinear2", trControl = trcontrol_path_SVM.feat,
                      metric = "Accuracy", data = pathCV_matrix)

varImp(fit_path_SVM.feat, scale = F)

# AUC
SVM.path.roc.feat <- evalm(fit_path_SVM.feat)
SVM.path.roc.feat$stdres


###################
## Decision tree ##
###################

set.seed(123)
trcontrol_path_dectree.feat <- trainControl(method = "LOOCV", savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
fit_path_dectree.feat <- train(condition2 ~ PWY.7237..myo...chiro..and.scyllo.inositol.degradation + PWY0.1221..putrescine.degradation.II +
                                 +RUMP.PWY..formaldehyde.oxidation.I, method = "C5.0", trControl = trcontrol_path_dectree.feat,
                          metric = "Accuracy", data = pathCV_matrix)

varImp(fit_path_dectree.feat, scale = F)

# AUC
dec.path.roc.feat <- evalm(fit_path_dectree.feat)
dec.path.roc.feat$stdres


####################
## Random forest ##
###################
set.seed(123)
trcontrol_path_rf.feat <- trainControl(method = "LOOCV", savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
fit_path_rf2.feat <- train(condition2 ~ PWY.7237..myo...chiro..and.scyllo.inositol.degradation + PWY0.1221..putrescine.degradation.II +
                             +RUMP.PWY..formaldehyde.oxidation.I, method = "rf", trControl = trcontrol_path_rf.feat,
                      metric = "Accuracy", data = pathCV_matrix)

varImp(fit_path_rf2.feat, scale = F)

# AUC
rf.path.roc.feat <- evalm(fit_path_rf2.feat)
rf.path.roc.feat$stdres






# gene family

##########
## k-NN ##
##########

set.seed(123)

trcontrol_gene_knn1.feat <- trainControl(method = "LOOCV",  savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
fit_gene_knn1.feat <- train(condition2 ~ UniRef90_A0A078ST01 + UniRef90_A0A139L4S8 +
                UniRef90_E2N7L9 , method = "knn", trControl = trcontrol_gene_knn1.feat,
                       metric = "Accuracy", data = geneCV_matrix)

varImp(fit_gene_knn1.feat, scale = F)

# AUC
knn.gene.roc.feat <- evalm(fit_gene_knn1.feat)
knn.gene.roc.feat$stdres


#################
## Naive Bayes ##
#################

trcontrol_gene_nb.feat <- trainControl(method = "LOOCV", savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
fit_gene_nb.feat <- train(condition2 ~ UniRef90_A0A078ST01 + UniRef90_A0A139L4S8 +
                       UniRef90_E2N7L9, method = "naive_bayes", trControl = trcontrol_gene_nb.feat,
                     metric = "Accuracy", data = geneCV_matrix)

varImp(fit_gene_knn1.feat, scale = F)

# AUC
nb.gene.roc.feat <- evalm(fit_gene_nb.feat)
nb.gene.roc.feat$stdres


################################
## Support vector machine SVM ##
################################
set.seed(123)
trcontrol_gene_SVM.feat <- trainControl(method = "LOOCV",savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
fit_gene_SVM.feat <- train(condition2 ~ UniRef90_A0A078ST01 + UniRef90_A0A139L4S8 +
                             UniRef90_E2N7L9, method = "svmLinear2", trControl = trcontrol_gene_SVM.feat,
                      metric = "Accuracy", data = geneCV_matrix)

varImp(fit_gene_SVM.feat, scale = F)

# AUC
SVM.gene.roc.feat <- evalm(fit_gene_SVM.feat)
SVM.gene.roc.feat$stdres

###################
## Decision tree ##
###################

set.seed(123)
trcontrol_gene_dectree.feat <- trainControl(method = "LOOCV", savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
fit_gene_dectree.feat <- train(condition2 ~ UniRef90_A0A078ST01 + UniRef90_A0A139L4S8 +
                            UniRef90_E2N7L9, method = "C5.0", trControl = trcontrol_gene_dectree.feat,
                          metric = "Accuracy", data = geneCV_matrix)

varImp(fit_gene_dectree.feat, scale = F)

# AUC
dec.gene.roc.feat <- evalm(fit_gene_dectree.feat)
dec.gene.roc.feat$stdres


####################
## Random forest ##
###################
set.seed(123)
trcontrol_gene_rf.feat <- trainControl(method = "LOOCV", savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
fit_gene_rf2.feat <- train(condition2 ~ UniRef90_A0A078ST01 + UniRef90_A0A139L4S8 +
                             UniRef90_E2N7L9, method = "rf", trControl = trcontrol_gene_rf.feat,
                      metric = "Accuracy", data = geneCV_matrix)

varImp(fit_gene_rf2.feat, scale = F)

# AUC
rf.gene.roc.feat <- evalm(fit_gene_rf2.feat)
rf.gene.roc.feat$stdres




