# Thema 1
library(GEOquery)
gse2553 <- getGEO('GSE2553',GSEMatrix=TRUE)
# show(gse2553)
# get Data & types of diseases (classes)
expressionData <- gse2553$GSE2553_series_matrix.txt.gz@assayData$exprs 
classes <- gse2553$GSE2553_series_matrix.txt.gz@phenoData@data$description

# create a dataframe having above data
# select only samples having below types
# 66 samples are returned 
result <-as.data.frame(cbind(classes, t(expressionData)))
in_data <-result[result$classes %in% c("Diagnosis: Liposarcoma", "Diagnosis: Malignant Fibrous Histiocytoma"), ]
class(in_data)

# select only the samples skip first column 
# these will be the input data

input_data <- dplyr::select(in_data, -1)
# use as a data frame
input_data <- as.data.frame(t(input_data))

# use below to display rows columns that have NA value
# which(is.na(input_data), arr.ind=TRUE)
# create a matrix with rows having a NA value somewhere
# test <- input_data[unique(which(is.na(input_data), arr.ind=TRUE)[,1]),]

# omit rows from input data having NA values
# get only numeric
all_input_data <- input_data
input_data <- na.omit(input_data)
input_data <- as.matrix(sapply(input_data, as.numeric))
# input data is 8459 x 66 
input_data_rows <- nrow(input_data)
input_data_cols <- ncol(input_data)
print(paste("genes found after removing NA values :", input_data_rows))
class(input_data)


# apply Rtsne to the transposed of the input data
# input data
library(Rtsne)
transp_input_data <- t(input_data)
tsne_out <- Rtsne(transp_input_data, perplexity = 20, dim=2) # Run TSNE
tsne_out$Y

# apply correlation to the transposed data
# cor_data matrix is a 8594x8594 matrix
cor_data <- cor(transp_input_data)

# save corr data into a tmp triangular matrix  set lower part and diagonal to 0
tmp <- cor_data
tmp[lower.tri(tmp)] <- 0
diag(tmp) <- 0

# exclude columns where correllation > 0.5
selected_data <- transp_input_data[, !apply(tmp, 2, function(x) any(abs(x) > 0.5, na.rm = TRUE))]
selected_data_rows <- nrow(selected_data)
selected_data_cols <- ncol(selected_data)
print(paste("genes found with corr < 0.5 :", selected_data_cols))
#head(selected_data)
tsne_out_corr <- Rtsne(selected_data, perplexity = 20, dim=2) # Run TSNE
tsne_out_corr$Y

# Thema 2
set.seed(0) 
library(factoextra)
library(cluster)
scaled_data <- scale(transp_input_data)
fviz_nbclust(scaled_data, kmeans, method = "wss")
