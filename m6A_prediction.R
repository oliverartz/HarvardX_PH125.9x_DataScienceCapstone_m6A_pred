## ----Setup, include = FALSE-----------------------------------------------------------------------------------
# caching chunks, so they don't run again when knitting report
knitr::opts_chunk$set(cache = TRUE, 
                      warning = FALSE, 
                      message = FALSE,
                      fig.align="center", 
                      out.width="70%",
                      library(dplyr))


## ----Load libraries, message=FALSE, warning=FALSE-------------------------------------------------------------
if (!require("tidyverse", quietly = TRUE)) 
  install.packages('tidyverse')
library(tidyverse)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
library(BiocManager)

if (!require("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) 
  install.packages("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)

if (!require("ChIPseeker", quietly = TRUE)) 
  install.packages("ChIPseeker")
library(ChIPseeker)

if (!require("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE)) 
  install.packages("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

if (!require("bedr", quietly = TRUE)) 
  install.packages('bedr')
library(bedr)

if (!require("caret", quietly = TRUE)) 
  install.packages('caret')
library(caret)

if (!require("doParallel", quietly = TRUE)) 
  install.packages('doParallel')
library(doParallel)

if (!require("bookdown", quietly = TRUE)) 
  install.packages('bookdown')
library(bookdown)

if (!require("formatR", quietly = TRUE)) 
  install.packages('formatR')
library(formatR)


## ----Load data, message=FALSE, warning=FALSE------------------------------------------------------------------
# fetching m6A site data from github repo
url <- "https://raw.githubusercontent.com/oliverartz/HarvardX_PH125.9x_DataScienceCapstone_m6A_pred/main/DARTseq_SuppTab2_Meyer_2019.csv"
df <- read_csv(url)

# cleanup
rm(url)


## ----Add non-methylated samples, message=FALSE, warning=FALSE-------------------------------------------------
# add column to specify methylated adenines
data <- df %>% dplyr::select(-`U/C`)
data$meth <- "methylated"

# generate approximately as many random adenines as in the data set
n_regions <- nrow(data) * 4

set.seed(2021, sample.kind = "Rounding")
random_regions <- get.random.regions(n = n_regions,
                   species = "human",
                   build = "hg19",
                   size.mean = 1,
                   size.sd = 0)

# add random strand information (with same distribution as in data set)
prob_plus <- table(data$Strand)[2] / (table(data$Strand)[1]+table(data$Strand)[2])
random_regions$strand <- sample(c("+","-"), 
                                n_regions, 
                                prob = c(prob_plus, 1-prob_plus), 
                                replace = TRUE)
random_regions$meth <- "not_methylated"

# filter for adenines
random_regions <- random_regions %>% 
  mutate(base = getSeq(Hsapiens, 
                       chr, 
                       start = start, 
                       end = start, 
                       strand = strand, 
                       as.character = TRUE)) %>% 
  filter(base == "A") %>% 
  dplyr::select(-base)

# add non-methylated regions to data frame
colnames(random_regions) <- colnames(data)
data <- rbind(data, random_regions)

# check for duplicates 
# to avoid labeling methylated adenines as unmethylated from the random set of loci
dupl <- data %>% 
  group_by(`C2U start`,Chr) %>% 
  summarize(n = n()) %>% 
  filter(n>1) %>% 
  pull(`C2U start`)

# remove duplicate
data <- data %>% filter(`C2U start` != dupl)

#cleanup
rm(random_regions, prob_plus,n_regions, dupl)


## ----Add features---------------------------------------------------------------------------------------------
# get genomic sequence adjacent to methylated adenine
data <- data %>% 
  mutate(seq_plus1 = getSeq(Hsapiens, 
                            Chr, 
                            start = `C2U start` + 1, 
                            end = `C2U start` + 1, 
                            strand = Strand, 
                            as.character = TRUE),
         seq_plus2 = getSeq(Hsapiens, 
                            Chr, start = `C2U start` + 2, 
                            end = `C2U start` + 2, 
                            strand = Strand, 
                            as.character = TRUE),
         seq_plus3 = getSeq(Hsapiens, 
                            Chr, 
                            start = `C2U start` + 3, 
                            end = `C2U start` + 3, 
                            strand = Strand, 
                            as.character = TRUE),
         seq_plus4 = getSeq(Hsapiens, 
                            Chr, 
                            start = `C2U start` + 4, 
                            end = `C2U start` + 4, 
                            strand = Strand, 
                            as.character = TRUE),
         seq_plus5 = getSeq(Hsapiens, 
                            Chr, 
                            start = `C2U start` + 5, 
                            end = `C2U start` + 5, 
                            strand = Strand, 
                            as.character = TRUE),
         seq_minus1 = getSeq(Hsapiens, 
                             Chr, 
                             start = `C2U start` - 1, 
                             end = `C2U start` - 1, 
                             strand = Strand, 
                             as.character = TRUE),
         seq_minus2 = getSeq(Hsapiens, 
                             Chr, 
                             start = `C2U start` - 2, 
                             end = `C2U start` - 2, 
                             strand = Strand, 
                             as.character = TRUE),
         seq_minus3 = getSeq(Hsapiens, 
                             Chr, 
                             start = `C2U start` - 3, 
                             end = `C2U start` - 3, 
                             strand = Strand, 
                             as.character = TRUE),
         seq_minus4 = getSeq(Hsapiens, 
                             Chr, 
                             start = `C2U start` - 4, 
                             end = `C2U start` - 4, 
                             strand = Strand, 
                             as.character = TRUE),
         seq_minus5 = getSeq(Hsapiens, 
                             Chr, 
                             start = `C2U start` - 5, 
                             end = `C2U start` - 5, 
                             strand = Strand, 
                             as.character = TRUE))

# get GC content of adjacent genomic sequence (proximal and distal)
proximal <- 50
distal <- 200

data <- data %>% 
  mutate(proximal_GC_3prime = letterFrequency(getSeq(Hsapiens, 
                                                     Chr, 
                                                     start = `C2U start` + 1, 
                                                     end = `C2U start` + proximal, 
                                                     strand = Strand), 
                                              letters = "GC", 
                                              as.prob = TRUE),
         proximal_GC_5prime = letterFrequency(getSeq(Hsapiens, 
                                                     Chr, 
                                                     start = `C2U start` - proximal, 
                                                     end = `C2U start` - 1, 
                                                     strand = Strand), 
                                              letters = "GC", 
                                              as.prob = TRUE),
         distal_GC_3prime = letterFrequency(getSeq(Hsapiens, 
                                                   Chr, 
                                                   start = `C2U start` + 1, 
                                                   end = `C2U start` + distal, 
                                                   strand = Strand), 
                                            letters = "GC", 
                                            as.prob = TRUE),
         distal_GC_5prime = letterFrequency(getSeq(Hsapiens, 
                                                   Chr, start = `C2U start` - distal, 
                                                   end = `C2U start` - 1, 
                                                   strand = Strand), 
                                            letters = "GC", 
                                            as.prob = TRUE))

# get feature annotation for m6A sites
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

peaks <- data %>% dplyr::select(1:4)
colnames(peaks) <- c("chr","start","end","strand")
peaks.gr <- makeGRangesFromDataFrame(peaks)

peak_annotation <- annotatePeak(peaks.gr, TxDb = txdb, verbose = FALSE)

peak_annotation <- peak_annotation %>% 
  as.data.frame() %>% 
  dplyr::select("annotation")

# clean up annotation (making it more general)
peak_annotation$annotation <- gsub("^Exon.*","Exon", peak_annotation$annotation)
peak_annotation$annotation <- gsub("^Promoter.*","Promoter", peak_annotation$annotation)
peak_annotation$annotation <- gsub("^Intron.*","Intron", peak_annotation$annotation)

# merge peak annotations with data frame
data <- cbind(data, peak_annotation)

# remove genomic locus from data frame
data <- data %>%  dplyr::select(-c(1:4))

# encode non-numeric columns as factors (necessary for some ML algorithms)
columns_to_factor <- c("seq_plus1", 
                       "seq_plus2", 
                       "seq_plus3", 
                       "seq_plus4", 
                       "seq_plus5", 
                       "seq_minus1", 
                       "seq_minus2", 
                       "seq_minus3", 
                       "seq_minus4", 
                       "seq_minus5", 
                       "annotation", 
                       "meth")

data[columns_to_factor] <- lapply(data[columns_to_factor], factor)

# cleanup 
rm(proximal, distal, txdb, peaks, peaks.gr, peak_annotation, columns_to_factor, ChIPseekerEnv)


## ----Data exploration - summary-------------------------------------------------------------------------------
# dimensions of data frame
dim <- dim(data)

# check for NAs
sum_NA <- sum(is.na(data))

# summary of data
summary_df <- summary(data)


## ----Initialize data frame model performance------------------------------------------------------------------
# initialize data frame for model performance
results_model_performance <- data.frame(model = character(),
                                        accuracy = numeric(),
                                        F1 = numeric(),
                                        wall_time_min = numeric())


## ----Split data set-------------------------------------------------------------------------------------------
set.seed(2022, sample.kind = "Rounding")
test_index <- createDataPartition(data$meth, times = 1, p = 0.2, list = FALSE)
test_set <- data %>% dplyr::slice(test_index)
train_set <- data %>% dplyr::slice(-test_index)

# cleanup
rm(test_index)


## ----Explore training and test test---------------------------------------------------------------------------
#proportion of methylated adenines in training
mean_a_train <- (mean(train_set$meth == "methylated") * 100) %>% round(digits = 3)

#proportion of methylated adenines in test
mean_a_test <- (mean(test_set$meth == "methylated") * 100) %>% round(digits = 3)


## ----Define cross validation parameters-----------------------------------------------------------------------
# cross validation parameters
fitControl <- trainControl(method = "repeatedcv",
                           number = 5, # n-fold cross-validation
                           repeats = 5, # repeat each cross-validation n times
                           classProbs = TRUE, # compute class probabilities
                           returnResamp = "final", # save resampled summary metrics
                           savePredictions = "final") # save the final predictions


## ----Set up multithreading------------------------------------------------------------------------------------
# define number of cores to use
cl <- makePSOCKcluster(5)

# register parallel backend
registerDoParallel(cl)


## ----Logistic regression--------------------------------------------------------------------------------------
# logistic regression model

# set seed to make results reproducible
set.seed(2022, sample.kind = "Rounding")

# model fitting
t_start <- Sys.time()
fit_glm <- train(meth ~ ., 
                 method = "glm", 
                 data = train_set,
                 trControl = fitControl)

t_end <- Sys.time()

# model validation
y_hat_glm <- predict(fit_glm, test_set)

# add results to performance data frame
cm <- confusionMatrix(y_hat_glm, mode = "everything", reference = test_set$meth)

results_model_performance[1, ] <- c(model = "logistic regression", 
                                    accuracy = cm$overall[1] %>% round(digits = 3),
                                    F1 = cm$byClass[7] %>% round(digits = 3),
                                    wall_time_min = difftime(t_end, t_start, units='min') %>% 
                                      as.numeric() %>% 
                                      round(digits = 2))

# cleanup
rm(cm, t_start, t_end)


## ----k-nearest neigbor model----------------------------------------------------------------------------------
# k-nearest neighbor model

# set seed to make results reproducible
set.seed(2022, sample.kind = "Rounding")

# define ks for optimization
k = data.frame(k = seq(2,20,2))

# model fitting
t_start <- Sys.time()
fit_knn <- train(meth ~ ., 
                 method = "knn", 
                 data = train_set, 
                 tuneGrid = k,
                 trControl = fitControl)
t_end <- Sys.time()

# validate model
y_hat_knn <- predict(fit_knn, test_set)

# add results to performance data frame
cm <- confusionMatrix(y_hat_knn, mode = "everything", reference = test_set$meth)

results_model_performance[2, ] <- c(model = "k-nearest neighbor", 
                                    accuracy = cm$overall[1] %>% round(digits = 3),
                                    F1 = cm$byClass[7] %>% round(digits = 3),
                                    wall_time_min = difftime(t_end, t_start, units='min') %>% 
                                      as.numeric() %>% 
                                      round(digits = 2))

# cleanup
rm(cm, t_start, t_end, k)


## ----Random forest model, message=FALSE, warning=FALSE--------------------------------------------------------
# random forest model with mtry optimization

# set seed to make results reproducible
set.seed(2022, sample.kind = "Rounding")

# use random search to find best mtry for random forest
t_start <- Sys.time()
# reduce number of trees to be produced to speed up optimization
ntree <- 5

# start fitting
rf_random <- train(meth ~ .,
                   data = train_set,
                   method = 'rf',
                   metric = 'Accuracy',
                   tuneLength  = 15, # number of randomly generated mtry values
                   trControl = fitControl,
                   ntree = ntree)

# pull best mtry value for subsequent model fitting
mtry <- data.frame(mtry = rf_random$bestTune %>% pull(mtry))

# set seed to make results reproducible
set.seed(2022, sample.kind = "Rounding")

# model fitting
fit_rf <- train(meth ~ ., 
                method = "rf", 
                data = train_set, 
                tuneGrid = mtry, 
                importance = TRUE,
                trControl = fitControl)
t_end <- Sys.time()

# validate model
y_hat_rf <- predict(fit_rf, test_set)

# add results to performance data frame
cm <- confusionMatrix(y_hat_rf, mode = "everything", reference = test_set$meth)

results_model_performance[3, ] <- c(model = "random forest", 
                                    accuracy = cm$overall[1] %>% round(digits = 3),
                                    F1 = cm$byClass[7] %>% round(digits = 3),
                                    wall_time_min = difftime(t_end, t_start, units='min') %>% 
                                      as.numeric() %>% 
                                      round(digits = 2))

# cleanup
rm(mtry, cm, t_start, t_end, ntree)


## ----Linear discriminant analysis-----------------------------------------------------------------------------
# linear discriminant analysis (LDA)

# set seed to make results reproducible
set.seed(2022, sample.kind = "Rounding")

# model fitting
t_start <- Sys.time()
fit_lda <- train(meth ~ ., 
                 method = "lda", 
                 data = train_set,
                 trControl = fitControl)
t_end <- Sys.time()

# model validation
y_hat_lda <- predict(fit_lda, test_set)

# add results to performance data frame
cm <- confusionMatrix(y_hat_lda, mode = "everything", reference = test_set$meth)

results_model_performance[4, ] <- c(model = "linear discriminant analysis", 
                                    accuracy = cm$overall[1] %>% round(digits = 3),
                                    F1 = cm$byClass[7] %>% round(digits = 3),
                                    wall_time_min = difftime(t_end, t_start, units='min') %>% 
                                      as.numeric() %>% 
                                      round(digits = 2))


# cleanup
rm(cm, t_start, t_end)


## ----Neural network, message=FALSE, warning=FALSE-------------------------------------------------------------
# neural network model

# set seed to make results reproducible
set.seed(2022, sample.kind = "Rounding")

# model fitting
t_start <- Sys.time()
fit_nn <- train(meth ~ ., 
                method = "nnet",
                data = train_set,
                trControl = fitControl,
                trace = FALSE)
t_end <- Sys.time()

# model validation
y_hat_nn <- predict(fit_nn, test_set)

# add results to performance data frame
cm <- confusionMatrix(y_hat_nn, mode = "everything", reference = test_set$meth)

results_model_performance[5, ] <- c(model = "neural network", 
                                    accuracy = cm$overall[1] %>% round(digits = 3),
                                    F1 = cm$byClass[7] %>% round(digits = 3),
                                    wall_time_min = difftime(t_end, t_start, units='min') %>% 
                                      as.numeric() %>% 
                                      round(digits = 2))

# cleanup
rm(cm, t_start, t_end)


## ----Make ensembl prediction----------------------------------------------------------------------------------
#accuracy of the ensemble model
prediction_matrix <- cbind(y_hat_glm, y_hat_knn, y_hat_lda, y_hat_rf,y_hat_nn)

votes <- rowMeans(prediction_matrix == 1)
y_hat_vote <- ifelse(votes > 0.5, 1, 0)
ensemble_accuracy <-(1 - mean(y_hat_vote == as.numeric(test_set$meth)-1)) %>% 
  round(digits = 3)

# re-factor votes from numeric to "methylated"/"not_methylated"
y_hat_vote <- y_hat_vote %>% 
  as.factor() %>% 
  recode_factor("1" = "methylated", "0" = "not_methylated")

# add results to performance data frame
cm <- confusionMatrix(y_hat_vote, mode = "everything", reference = test_set$meth)

results_model_performance[6, ] <- c(model = "ensemble", 
                                    accuracy = cm$overall[1] %>% round(digits = 3),
                                    F1 = cm$byClass[7] %>% round(digits = 3),
                                    wall_time_min = NA)


#which is the best model?
accuracy <- colMeans(prediction_matrix == as.numeric(test_set$meth))


# cleanup
rm(accuracy, prediction_matrix, votes, y_hat_vote, ensemble_accuracy)


## ----Preparing Linder (2015) data-----------------------------------------------------------------------------
# fetching m6A site data from github repo
url <- "https://raw.githubusercontent.com/oliverartz/HarvardX_PH125.9x_DataScienceCapstone_m6A_pred/main/CITS_miCLIP_SuppTab1_Linder_2015.csv"
df2 <- read_csv(url)

# cleanup
rm(url)

# select relevant columns
data2 <- df2 %>% dplyr::select(2,3,4,6,5)

# copy column names from other data frame to make later steps easier
colnames(data2) <- colnames(df)

# add column to specify methylated adenines
data2$meth <- "methylated"

# generate approximately as many random adenines as in the data set
n_regions <- nrow(data2) * 4

set.seed(2021, sample.kind = "Rounding")
random_regions <- get.random.regions(n = n_regions,
                   species = "human",
                   build = "hg19",
                   size.mean = 1,
                   size.sd = 0,
                   )

# add random strand information (with same distribution as in data set)
prob_plus <- table(data2$Strand)[2] / (table(data2$Strand)[1]+table(data2$Strand)[2])
random_regions$strand <- sample(c("+","-"), 
                                n_regions, 
                                prob = c(prob_plus, 1-prob_plus), 
                                replace = TRUE)
random_regions$meth <- "not_methylated"

# filter for adenines
random_regions <- random_regions %>% 
  mutate(base = getSeq(Hsapiens, 
                       chr, 
                       start = start, 
                       end = start, 
                       strand = strand, 
                       as.character = TRUE)) %>% 
  filter(base == "A") %>% 
  dplyr::select(-base)

random_regions <- random_regions %>% 
  mutate(fill = 1) %>% 
  dplyr::select(1,2,3,6,4,5)

# add non-methylated regions to data frame
colnames(random_regions) <- colnames(data2)
data2 <- rbind(data2, random_regions)

# check for duplicates 
# to avoid labeling methylated adenines as unmethylated from the random set of loci
dupl <- data2 %>% 
  group_by(`C2U start`,Chr) %>% 
  summarize(n = n()) %>% 
  filter(n>1) %>% 
  pull(`C2U start`)

# remove duplicate
#data2 <- data2 %>% filter(`C2U start` != dupl)

#cleanup
rm(random_regions, prob_plus,n_regions, dupl)

# get genomic sequence adjacent to methylated adenine
data2 <- data2 %>% 
  mutate(seq_plus1 = getSeq(Hsapiens, 
                            Chr, 
                            start = `C2U start` + 1, 
                            end = `C2U start` + 1, 
                            strand = Strand, 
                            as.character = TRUE),
         seq_plus2 = getSeq(Hsapiens, 
                            Chr, 
                            start = `C2U start` + 2,
                            end = `C2U start` + 2, 
                            strand = Strand, 
                            as.character = TRUE),
         seq_plus3 = getSeq(Hsapiens, 
                            Chr, 
                            start = `C2U start` + 3, 
                            end = `C2U start` + 3, 
                            strand = Strand, 
                            as.character = TRUE),
         seq_plus4 = getSeq(Hsapiens, 
                            Chr, 
                            start = `C2U start` + 4, 
                            end = `C2U start` + 4, 
                            strand = Strand, 
                            as.character = TRUE),
         seq_plus5 = getSeq(Hsapiens, 
                            Chr, 
                            start = `C2U start` + 5, 
                            end = `C2U start` + 5, 
                            strand = Strand, 
                            as.character = TRUE),
         seq_minus1 = getSeq(Hsapiens, 
                             Chr, 
                             start = `C2U start` - 1, 
                             end = `C2U start` - 1, 
                             strand = Strand, 
                             as.character = TRUE),
         seq_minus2 = getSeq(Hsapiens, 
                             Chr, 
                             start = `C2U start` - 2, 
                             end = `C2U start` - 2, 
                             strand = Strand, 
                             as.character = TRUE),
         seq_minus3 = getSeq(Hsapiens, 
                             Chr, 
                             start = `C2U start` - 3, 
                             end = `C2U start` - 3, 
                             strand = Strand, 
                             as.character = TRUE),
         seq_minus4 = getSeq(Hsapiens, 
                             Chr, 
                             start = `C2U start` - 4, 
                             end = `C2U start` - 4, 
                             strand = Strand, 
                             as.character = TRUE),
         seq_minus5 = getSeq(Hsapiens, 
                             Chr, 
                             start = `C2U start` - 5, 
                             end = `C2U start` - 5, 
                             strand = Strand, 
                             as.character = TRUE))

# get GC content of adjacent genomic sequence (proximal and distal)
proximal <- 50
distal <- 200

data2 <- data2 %>% 
  mutate(proximal_GC_3prime = letterFrequency(getSeq(Hsapiens, 
                                                     Chr, 
                                                     start = `C2U start` + 1, 
                                                     end = `C2U start` + proximal, 
                                                     strand = Strand), 
                                              letters = "GC", 
                                              as.prob = TRUE),
         proximal_GC_5prime = letterFrequency(getSeq(Hsapiens, 
                                                     Chr, 
                                                     start = `C2U start` - proximal,
                                                     end = `C2U start` - 1, 
                                                     strand = Strand), 
                                              letters = "GC", 
                                              as.prob = TRUE),
         distal_GC_3prime = letterFrequency(getSeq(Hsapiens, 
                                                   Chr, 
                                                   start = `C2U start` + 1, 
                                                   end = `C2U start` + distal, 
                                                   strand = Strand), 
                                            letters = "GC", 
                                            as.prob = TRUE),
         distal_GC_5prime = letterFrequency(getSeq(Hsapiens,
                                                   Chr, 
                                                   start = `C2U start` - distal, 
                                                   end = `C2U start` - 1, 
                                                   strand = Strand), 
                                            letters = "GC", 
                                            as.prob = TRUE))

# get feature annotation for m6A sites
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

peaks <- data2 %>% dplyr::select(1,2,3,5)
colnames(peaks) <- c("chr","start","end","strand")
peaks.gr <- makeGRangesFromDataFrame(peaks)

peak_annotation <- annotatePeak(peaks.gr, TxDb = txdb, verbose = FALSE)

peak_annotation <- peak_annotation %>% 
  as.data.frame() %>% 
  dplyr::select(1,2,3,5,6)
colnames(peak_annotation) <- c("Chr", 
                               "C2U start", 
                               "C2U end", 
                               "Strand",
                               "annotation")

# clean up annotation (making it more general)
peak_annotation$annotation <- gsub("^Exon.*","Exon", peak_annotation$annotation)
peak_annotation$annotation <- gsub("^Promoter.*","Promoter", peak_annotation$annotation)
peak_annotation$annotation <- gsub("^Intron.*","Intron", peak_annotation$annotation)

# merge peak annotations with data frame
data2 <- left_join(data2, 
                   peak_annotation, 
                   by = c("Chr", "C2U start", "C2U end", "Strand"))

# annotating mitochondrial genes leads to incomplete cases
# filter incomplete cases / mitochondrial genes
data2 <- data2 %>% na.omit

# remove genomic locus from data frame
data2 <- data2 %>% dplyr::select(-c(1:5))

# encode non-numeric columns as factors (necessary for some ML algorithms)
columns_to_factor <- c("seq_plus1", 
                       "seq_plus2", 
                       "seq_plus3", 
                       "seq_plus4", 
                       "seq_plus5", 
                       "seq_minus1", 
                       "seq_minus2", 
                       "seq_minus3", 
                       "seq_minus4", 
                       "seq_minus5", 
                       "annotation", 
                       "meth")

data2[columns_to_factor] <- lapply(data2[columns_to_factor], factor)


# cleanup
rm(distal, proximal, txdb, peak_annotation, peaks.gr)


## ----Testing models on Linder (2015) data---------------------------------------------------------------------
# applying fits to new data
y_hat_glm_Linder <- predict(fit_glm, data2) 
y_hat_knn_Linder <- predict(fit_knn, data2) 
y_hat_lda_Linder <- predict(fit_lda, data2) 
y_hat_rf_Linder <- predict(fit_rf, data2)
y_hat_nn_Linder <- predict(fit_nn, data2)

# accuracy of the ensemble model
prediction_matrix <- cbind(y_hat_glm_Linder, 
                           y_hat_knn_Linder, 
                           y_hat_rf_Linder, 
                           y_hat_lda_Linder,
                           y_hat_nn_Linder)

votes <- rowMeans(prediction_matrix == 1)
y_hat_vote <- ifelse(votes > 0.5, 1, 0)

ensemble_accuracy <- (1 - mean(y_hat_vote == as.numeric(data2$meth)-1)) %>% 
  round(digits = 3)

# re-factor votes from numeric to "methylated"/"not_methylated"
y_hat_vote <- y_hat_vote %>% 
  as.factor() %>% 
  recode_factor("1" = "methylated", "0" = "not_methylated")


# add accuracy on Linder to results data frame
cm <- confusionMatrix(y_hat_vote, mode = "everything", reference = data2$meth)

results_model_performance_Linder <- data_frame(model = "ensemble", 
                                    accuracy = cm$overall[1] %>% round(digits = 3),
                                    F1 = cm$byClass[7] %>% round(digits = 3))
# cleanup
rm(accuracy_Lin, prediction_matrix, votes, y_hat_vote, ensemble_accuracy)


## ----tab_performance_meyer------------------------------------------------------------------------------------
results_model_performance %>% 
  rename("model" = "model", 
         "accuracy" = "accuracy", 
         "F1" = "F1", 
         "wall_time_min" = "wall time (min)") %>% 
  knitr::kable(caption = "Model performance on Meyer data",
               align = "lrrr")


## ----parameter tuning knn, fig.cap="Tuning number of neighbors (k) for knn model"-----------------------------
plot(fit_knn)


## ----parameter tuning mtry, fig.cap="Tuning number of considered factors (mtry) for Random Forest model"------
plot(rf_random)


## ----Assess variable importance, out.width = "48%", fig.show = "hold", fig.ncol = 2, fig.cap = "Variable importance"----
# number of variables to plot
n_vars <- 10

# make plots
plot(varImp(fit_glm), 
     top = n_vars, 
     adj = 0, 
     main = "(A) Linear regression")

plot(varImp(fit_knn), 
     top = n_vars,
     adj = 0,
     main = "(B) K-nearest neigbors")

plot(varImp(fit_nn), 
     top = n_vars, 
     adj = 0, 
     main = "(C) Neural network")

plot(varImp(fit_rf), 
     top = n_vars, 
     adj = 0, 
     main = "(D) Random forest")

plot(varImp(fit_lda), 
     top = n_vars, 
     adj = 0, 
     main = "(E) Linear discriminant analysis")

# cleanup
rm(n_vars)


## ----Session info---------------------------------------------------------------------------------------------
sessionInfo()

