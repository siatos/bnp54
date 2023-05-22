# Thema 1
library(GEOquery)
gse2553 <- getGEO('GSE2553',GSEMatrix=TRUE)
# show(gse2553)
# get Data & types of diseases (classes)
# expressionData is 12600(genes) x 181(samples)
expressionData <- gse2553$GSE2553_series_matrix.txt.gz@assayData$exprs 
classes <- gse2553$GSE2553_series_matrix.txt.gz@phenoData@data$description

# create a dataframe having above data [in a matrix where samples are the rows and columns are the genes (characteristics/features)]
# select only samples having below types
# 66 samples are returned : a matrix 66x12601
result <-as.data.frame(cbind(classes, t(expressionData)))
in_data <-result[result$classes %in% c("Diagnosis: Liposarcoma", "Diagnosis: Malignant Fibrous Histiocytoma"), ]
class(in_data)

# select only the samples skip first column 
# these will be the input data

input_data <- dplyr::select(in_data, -1)
# use as a data frame
# get transpose 12600x66 
input_data <- as.data.frame(t(input_data))

# use below to display rows columns that have NA value
# which(is.na(input_data), arr.ind=TRUE)
# create a matrix with rows having a NA value somewhere
# test <- input_data[unique(which(is.na(input_data), arr.ind=TRUE)[,1]),]

# omit rows from input data having NA values
# get only numeric
input_data <- na.omit(input_data)
input_data <- as.matrix(sapply(input_data, as.numeric))
# now input data is 8459 (genes) x 66 (samples) 
input_data_rows <- nrow(input_data)
input_data_cols <- ncol(input_data)
print(paste("genes found after removing NA values :", input_data_rows))
class(input_data)


# apply Rtsne to the transposed of the input data so that samples (66) are the rows
# 66x8594
library(Rtsne)
transp_input_data <- t(input_data)
tsne_out <- Rtsne(transp_input_data, perplexity = 20, dim=2) # Run TSNE
tsne_out$Y
plot(tsne_out$Y[, 1] ~ tsne_out$Y[, 2], pch = 10, col = "black")

# apply correlation to the transposed data
# cor_data matrix is a 8594x8594 matrix
cor_data <- cor(transp_input_data)

# save corr data into a tmp triangular matrix  set upper part and diagonal to 0
tmp <- cor_data
tmp[upper.tri(tmp)] <- 0
diag(tmp) <- 0

# exclude columns where correllation > 0.5
# selected_data <- transp_input_data[, !apply(tmp, 2, function(x) any(abs(x) > 0.5, na.rm = TRUE))]
selected_data <- transp_input_data[, !apply(tmp, 2, function(x) any(x > 0.5, na.rm = TRUE))]
selected_data_rows <- nrow(selected_data)
selected_data_cols <- ncol(selected_data)
# selected data is 66x462
print(paste("applying corr<0.5 "))
print(paste("genes found with corr < 0.5 :", selected_data_cols))

# apply Rtsne to the selected of the input data (corr <= 0.5)
tsne_out_corr <- Rtsne(selected_data, perplexity = 20, dim=2) # Run TSNE
tsne_out_corr$Y
plot(tsne_out$Y[, 1] ~ tsne_out$Y[, 2], pch = 10, col = "red")


# Thema 2
set.seed(1) 
library(factoextra)
library(cluster)

# scaled selected data (66x462 scaled)
# (66 samples x 462 genes)
scaled_data <- scale(transp_input_data)
class(scaled_data)
scaled_data_rows <- nrow(scaled_data)
scaled_data_cols <- ncol(scaled_data)
fviz_nbclust(scaled_data, kmeans, method = "wss", k.max = 15)
km_selected_scaled <- kmeans(scaled_data, centers = 11, nstart = 25)
fviz_cluster(km_selected_scaled, data = scaled_data)


# non-scaled selected data (66x462 non scaled)
fviz_nbclust(transp_input_data, kmeans, method = "wss", k.max = 15)
km_selected_non_scaled <- kmeans(transp_input_data, centers = 11, nstart = 25)
fviz_cluster(km_selected_non_scaled, data = transp_input_data)


# original scaled data
origin_data <- as.data.frame(result[, !names(result) %in% c("classes")]) 
origin_data1 <- na.omit(t(origin_data))
origin_data2 <- as.matrix(sapply(origin_data1, as.numeric))
scaled_origin_data <- scale(origin_data2)
fviz_nbclust(scaled_origin_data, kmeans, method = "wss", k.max = 10)
km_scaled_origin_data <- kmeans(scaled_origin_data, centers = 4, nstart = 25)
fviz_cluster(km_scaled_origin_data, data = scaled_origin_data)


#head(scaled_data)

#te <- dplyr::select(result, -1)
#te <- as.data.frame(t(te))
#te <- na.omit(te)
#class(te)
#te <- as.matrix(sapply(te, as.numeric))
#scaled_data <- scale(te)

# fviz_nbclust(scaled_data, kmeans, method = "wss")  
#  geom_vline(xintercept = 3, linetype = 2)

# Thema 3
# need to install
BiocManager::install("randomForest")
