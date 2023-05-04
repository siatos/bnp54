# full path for the needed input file
# needs to be adjusted accordingly
exercise_txt_path <- "C:/Users/siatos/Desktop/pae/bnp54/GE/GE2/Rcode/GE2_3/exercise2.txt"

################## θέμα 3: Q1 a-c  ###############################################
library(readr)
#read file input data in a data frame
input_data_df <- read.table(exercise_txt_path, sep="\t", header=TRUE)
class(input_data_df)
## get gene names as a vector
gene_names <- as.matrix(input_data_df$ProbeID)

# data frame - keep only numerical data
input_data <- input_data_df[, unlist(lapply(input_data_df, is.numeric), use.names = FALSE)]
rownames(input_data) <- gene_names
class(input_data)

# get first 100 rows from input data
input_data_100 <- head(input_data, n=100)
class(input_data_100)
## convert data frame of first 100 into a numeric matrix and get heatmap
## (for any heatmap usually matrix = genes(rows) X samples(columns))
heatmap(as.matrix(sapply(input_data_100, as.numeric)))
boxplot(input_data_100, xlab="samples ", ylab="values in Genes", main="sample boxplot - only 100 genes")
#boxplot(input_data, xlab="samples ", ylab="values in Genes", main="sample boxplot - only 100 genes")

###################################################################################

################## θέμα 3: Q2 a-b  ################################################
## normalization - scale to normal distribution  Z(0,1)
## normalize using scale we need to scale per gene since scale works on columns 
## transpose input data: genes will become columns and then apply scale and re transpose (genes will be rows) for boxplot 
tinput_data <- t(input_data)
class(tinput_data)
norm_data <- as.data.frame(t(scale(tinput_data, center=TRUE, scale=TRUE)))
class(norm_data)
boxplot(norm_data, xlab="samples", ylab="values", main="boxplot z-norm data")

#normalization - log
log_data <- log2(input_data)


## applying log2 gives -Inf result if value is close  to 0 these rows need to be removed
## apply rowwise (MARGIN=1) function is.finite to get only rows with finite values
log_data <- as.data.frame(log_data[apply(log_data, 1, function(x) all(is.finite(x))), ])
class(log_data)
boxplot(log_data, xlab="samples", ylab="values", main="boxplot log data")

## normalization - scale in [0, 1]
## apply min=01/max=1 normalization on every row (gene)
## apply funcion transposes the original data
## so we need to re-transpose for boxplot (per sample)
class(input_data)
scale_data <- apply(input_data, 1, function(x)(x-min(x))/(max(x)-min(x)))
class(scale_data)
## scale_data is 8x12626 now, so for boxplot per sample we need 12626x8 
boxplot(as.data.frame(t(scale_data)), xlab="samples", ylab="values", main="boxplot scale data")
###################################################################################


################## θέμα 3: Q3 a-b  ################################################
log_data_adult <- log_data[, c("brain.1", "brain.2", "liver.1", "liver.2")]
log_data_fetal <- log_data[, c("fetal.brain.1", "fetal.brain.2", "fetal.liver.1", "fetal.liver.2")]

log_data_adult_mean <- apply(log_data_adult, 1, mean)
log_data_fetal_mean <- apply(log_data_fetal, 1, mean)
log_data_mean <- cbind(log_data_adult_mean, log_data_fetal_mean)
colnames(log_data_mean) <- c('adult', 'fetal')
log_data_mean <- data.frame(log_data_mean)
log2FC <- data.frame(log_data_mean$adult - log_data_mean$fetal)
colnames(log2FC) <- c('log2FC')


## extended log data: all columns plus means for each gene - two groups and respective log2FC value total 10 columns 
ext_log_data <- data.frame(cbind(log_data_adult, log_data_fetal, log_data_mean, abs(log2FC)))

## select all rows with abs(log2FC) > 2 and order in desc order: from 12616 -> 525 objects
ext_log_data <- dplyr::filter(ext_log_data, log2FC > 2) 
ext_log_data <- dplyr::arrange(ext_log_data, desc(log2FC))
class(ext_log_data)
colnames(ext_log_data)

## select all data columns except means and log2FC for heatmap only
heat_log_data <- ext_log_data[c("brain.1", "brain.2", "liver.1", "liver.2", "fetal.brain.1", "fetal.brain.2", "fetal.liver.1", "fetal.liver.2")]
heat_log_data <- as.matrix(heat_log_data)

## heatmap of (sorted in desc order) log2 data with abs(log2FC) > 2
heatmap(heat_log_data)

## t-test
ncols <- ncol(ext_log_data)
nrows <- nrow(ext_log_data)

pvalues <- matrix(nrow=nrows, ncol=1)
for(i in 1:nrows) {
  tt<- t.test(ext_log_data[i,1:4], ext_log_data[i, 5:8])
  pvalues[i, 1] <- tt$p.value
}

pvalues <- as.data.frame(p.adjust(pvalues, method = "fdr"))
colnames(pvalues) <- c('pvalues')

ext_log_data_with_adjusted <- data.frame(cbind(ext_log_data, pvalues))
ext_log_data_with_adjusted <- dplyr::filter(ext_log_data_with_adjusted, pvalues < 0.05)
ext_log_data_with_adjusted

## select all data columns for heatmap only
heat_log_data_adjusted <- as.matrix(ext_log_data_with_adjusted[c("brain.1", "brain.2", "liver.1", "liver.2", "fetal.brain.1", "fetal.brain.2", "fetal.liver.1", "fetal.liver.2")])

## heatmap of (sorted in desc order) log2 data with abs(log2FC) > 2 and fdr < 0.05
heatmap(heat_log_data_adjusted)

###################################################################################



################## θέμα 3: Q4 a-d  ################################################
## anova
new_log_data <- as.matrix(log_data[c("brain.1", "brain.2", "liver.1", "liver.2", "fetal.brain.1", "fetal.brain.2", "fetal.liver.1", "fetal.liver.2")])
anova_pval <- matrix(nrow=nrow(new_log_data), ncol=1) #Produce vector to save p-value

anova_data <- list()
for (i in 1:nrow(new_log_data)) { 
  data <- data.frame(group = rep(c(1, 1, 2, 2, 3, 3, 4, 4), each = 1),
                     values = new_log_data[i, ])
  summary_data <- summary(aov(values~factor(group), data=data))
  anova_data[[i]] <- summary_data
  anova_pval[i, 1] <- summary_data[[1]][["Pr(>F)"]][[1]]
}

rownames(anova_pval) <- rownames(new_log_data)
colnames(anova_pval) <- c('anova_pval')

## create a data frame with gene values and calculated p value there will 12616 rows
anova_vals <- data.frame(cbind(new_log_data, anova_pval))

## select only these with p < 000001
## and plot heatmap
anova_vals <- dplyr::filter(anova_vals, anova_pval < 0.000001)
anova_heatmap_vals <- as.matrix(anova_vals[c("brain.1", "brain.2", "liver.1", "liver.2", "fetal.brain.1", "fetal.brain.2", "fetal.liver.1", "fetal.liver.2")])
heatmap(anova_heatmap_vals)

## select only 5 genes p < 0.0000005
anova_5s <- dplyr::filter(anova_vals, anova_pval < 0.0000005)

## Tukey 
new_5_data <- as.matrix(anova_5s[c("brain.1", "brain.2", "liver.1", "liver.2", "fetal.brain.1", "fetal.brain.2", "fetal.liver.1", "fetal.liver.2")])

anova_5_data <- list()
for (i in 1:nrow(new_5_data)) { 
  data <- data.frame(group = rep(c(1, 1, 2, 2, 3, 3, 4, 4), each = 1),
                     values = new_5_data[i, ])
  ## get anova results for th5 selected genes
  model_data <- aov(values~factor(group), data=data)
  ## run TukeyHSD  
  res <- TukeyHSD(model_data, conf.level=.95)
  ## plot TukeyHSD
  plot(res, las = 2)
}

###################################################################################

################## θέμα 3: Q5 a-b  ################################################
## pca
# use all data (12626x8): input_data above 
# get transpose 8 samples X 12626 characteristics (genes)
input_data
pca_all_data <- t(as.matrix(input_data))
pca_results_all <- prcomp(pca_all_data, scale = TRUE)
summary(pca_results_all)

library(ggfortify)
autoplot(pca_results_all)

## use reduced data with abs(log2FC) > 2 : heat_log_data  above
## get transpose 8 samples X 525 characteristics (genes)
## 
pca_log_data <- t(as.matrix(heat_log_data))
pca_results_log_data <- prcomp(pca_log_data, scale = TRUE)
summary(pca_results_log_data)
library(ggfortify)
autoplot(pca_results_log_data)

###################################################################################






