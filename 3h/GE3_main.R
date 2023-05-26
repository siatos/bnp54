## Thema 1
library(GEOquery)
gse2553 <- getGEO('GSE2553',GSEMatrix=TRUE)
# show(gse2553)
## get Data & types of diseases (classes)
## expressionData is 12600(genes) x 181(samples)
expressionData <- gse2553$GSE2553_series_matrix.txt.gz@assayData$exprs 
classes <- gse2553$GSE2553_series_matrix.txt.gz@phenoData@data$description

## create a dataframe having above data [in a matrix where samples are the rows and columns are the genes (characteristics/features)]
## select only samples having below types
## 66 samples are returned : a matrix 66x12601
result <-as.data.frame(cbind(classes, t(expressionData)))
in_data <-result[result$classes %in% c("Diagnosis: Liposarcoma", "Diagnosis: Malignant Fibrous Histiocytoma"), ]
class(in_data)

# select only the samples skip first column 
# these will be the input data
input_data <- dplyr::select(in_data, -1)

## use as a data frame
## get transpose 12600x66 
input_data <- as.data.frame(t(input_data))

## note: use below to display rows columns that have NA value
## which(is.na(input_data), arr.ind=TRUE)
## create a matrix with rows having a NA value somewhere
## test <- input_data[unique(which(is.na(input_data), arr.ind=TRUE)[,1]),]

##
##  omit rows from input data having NA values
##  na.omit removes ROWs on columns having na values
##  so we need to have genes in the ROWs
##  get only numeric
## 
input_data <- na.omit(input_data)
input_data <- as.matrix(sapply(input_data, as.numeric))

## now input data is 8594 (genes) x 66 (samples) 
input_data_rows <- nrow(input_data)
input_data_cols <- ncol(input_data)
print(paste("genes found after removing NA values :", input_data_rows))
class(input_data)
str(input_data)


## apply Rtsne to the transposed of the input data so that samples (66) are the rows
## 66x8594
library(Rtsne)
transp_input_data <- t(input_data)
row_names <- rownames(transp_input_data)
row_names
tsne_out <- Rtsne(transp_input_data, perplexity = 20, dim=2) # Run TSNE
tsne_out$Y
plot(tsne_out$Y[, 1] ~ tsne_out$Y[, 2], pch = 20, col = "black")

## apply correlation to the transposed data
## cor_data matrix is a 8594x8594 matrix
cor_data <- cor(transp_input_data)
## save corr data into a tmp triangular matrix  set upper part and diagonal to 0
tmp <- cor_data
tmp[upper.tri(tmp)] <- 0
diag(tmp) <- 0

## exclude columns where correllation > 0.5
selected_data <- transp_input_data[, !apply(tmp, 2, function(x) any(abs(x) > 0.5, na.rm = TRUE))]
## selected_data <- transp_input_data[, !apply(tmp, 2, function(x) any(x > 0.5, na.rm = TRUE))]
selected_data_rows <- nrow(selected_data)
selected_data_cols <- ncol(selected_data)
## selected data is 66x413
print(paste("applying corr<0.5 ", "genes found with corr < 0.5 :", selected_data_cols))

## apply Rtsne to the selected of the input data (corr <= 0.5)
tsne_out_corr <- Rtsne(selected_data, perplexity = 20, dim=2) # Run TSNE
tsne_out_corr$Y
plot(tsne_out$Y[, 1] ~ tsne_out$Y[, 2], pch = 20, col = "red")


## Thema 2
set.seed(1) 
## First case: scale selected data (66x8459)
## 
scaled_data <- scale(transp_input_data)
class(scaled_data)

scaled_data_rows <- nrow(scaled_data)
scaled_data_cols <- ncol(scaled_data)
print(paste("scaled selected data, rows:", scaled_data_rows, " columns:", scaled_data_cols))

print("Find best kmeans cluster separation  k = 1..12 scaled selected data")
scaled_result <- factoextra::fviz_nbclust(scaled_data, kmeans, method = "wss", k.max = 12)
scaled_result
scaled_result$data

## run k-means for scale data & cluster centers
km_selected_scaled <- kmeans(scaled_data, centers = 7, nstart = 25)

## display clusters
clusters_inf_scaled <- factoextra::fviz_cluster(km_selected_scaled, data = scaled_data)
clusters_inf_scaled

## Second case: non-scaled selected data (66x8459)
##
print("Find best kmeans cluster separation k = 1..12 non-scaled selected data")
non_scaled_result <- factoextra::fviz_nbclust(transp_input_data, kmeans, method = "wss", k.max = 12)
non_scaled_result
non_scaled_result$data

## run k-means for non scaled data & cluster centers
km_selected_non_scaled <- kmeans(transp_input_data, centers = 7, nstart = 25)
clusters_inf_non_scaled <- factoextra::fviz_cluster(km_selected_non_scaled, data = transp_input_data)
clusters_inf_non_scaled


##
# Dissimilarity matrix
d <- dist(scaled_data, method = "euclidean")
# Hierarchical clustering using single Linkage
hc1 <- hclust(d, method = "single" )

# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1,  main="Dendrogram Scaled Data - samples 66x8594")

## run agnes: Hierarchical Agglomerative with single linkage
##res_agnes <- cluster::agnes(x = scaled_data, # data matrix
##                            stand = TRUE, # Standardize the data
##                            metric = "euclidean", # metric for distance matrix
##                            method = "single") # Linkage method

#plot(res_agnes, labels = cS, nmax = 150)# bannerplot labels are mess
#plot(res_agnes, cex = 0.6, hang = -1)
scaled_tsne <- as.data.frame(scale(tsne_out$Y))
rownames(scaled_tsne) <- row_names
factoextra::fviz_nbclust(scaled_tsne, kmeans, method = "wss", k.max = 12)
km_sne <- kmeans(scaled_tsne, centers = 5, nstart = 25)
factoextra::fviz_cluster(km_sne, data = scaled_tsne)

d <- dist(scaled_data, method = "euclidean")
hc1 <- hclust(d, method = "single" )

# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1, main="Dendrogram sne Data -samples 66 x 2")




# Thema 3
# need to install
#BiocManager::install("randomForest")

