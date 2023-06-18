## full path for the needed input file
## needs to be (can be) adjusted accordingly
## here we use "root" folder c:/

exercise_txt_path <- "C:/data_bnp54.txt"
library(readr)
## read file input data in a data frame
in_data <- read.table(exercise_txt_path, sep="\t", header=TRUE)
class(in_data)
in_data
str(in_data)
colnames <- colnames(in_data)

# Load libraries
library(cluster)
library(factoextra)

set.seed(4300) 
## scale selected data (502x100)
scaled_data <- sapply(in_data, as.numeric)
scaled_data <- scale(scaled_data)
rownames(scaled_data) <- row.names(in_data)
colnames(scaled_data) <- colnames
print(paste("scaled selected data, rows:", nrow(scaled_data), " columns:", ncol(scaled_data)))

## run kmeans for scaled data & cluster centers
## select k=2  
##
km_res2 <- kmeans(scaled_data, centers = 2, nstart = 25)
## display clusters
kmclusters2 <- factoextra::fviz_cluster(km_res2, data = scaled_data)
kmclusters2

## run kmeans for scaled data & cluster centers
## select k=4  
##
km_res4 <- kmeans(scaled_data, centers = 4, nstart = 25)
## display clusters
kmclusters4 <- factoextra::fviz_cluster(km_res4, data = scaled_data)
kmclusters4

tsne_data <- as.matrix(sapply(in_data, as.numeric))
colnames(tsne_data) <- colnames
library(Rtsne)
tsne_out <- Rtsne(tsne_data, perplexity = 20, dim=2) # Run TSNE
tsne_plot <- data.frame(x = tsne_out$Y[ ,1], y = tsne_out$Y[, 2]) 
library(ggplot2)
plot(tsne_out$Y, main="Gene Expression Data", sub="Reduced to 2 features (t-SNE)", xlab="t-SNE 1", ylab="t-SNE 2")

# tsne with kmeans
## run kmeans on tsne 1st dim
## select k=2  
##
tsne_km_res2 <- kmeans(tsne_out$Y[,1], centers = 2, nstart = 25)
## display clusters
tsne_kmclusters2 <- factoextra::fviz_cluster(tsne_km_res2, data = scaled_data)
tsne_kmclusters2

## run kmeans on tsne 1st dim
## select k=4  
##
tsne_km_res4 <- kmeans(tsne_out$Y[,1], centers = 4, nstart = 25)
## display clusters
tsne_kmclusters4 <- factoextra::fviz_cluster(tsne_km_res4, data = scaled_data)
tsne_kmclusters4



# Compute silhouette values
# 1
sil_km_res2 <- silhouette(km_res2$cluster, dist(scaled_data))
fviz_silhouette(sil_km_res2)

# 2
sil_km_res4 <- silhouette(km_res4$cluster, dist(scaled_data))
fviz_silhouette(sil_km_res4)

# 3
sil_tsne_res2 <- silhouette(tsne_km_res2$cluster, dist(tsne_out$Y[,1]))
fviz_silhouette(sil_tsne_res2)

# 4
sil_tsne_res4 <- silhouette(tsne_km_res4$cluster, dist(tsne_out$Y[,1]))
fviz_silhouette(sil_tsne_res4)

