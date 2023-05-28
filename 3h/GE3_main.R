## Thema 1
library(GEOquery)
gse2553 <- getGEO('GSE2553',GSEMatrix=TRUE)

## get Data & types of diseases (classes)
## expressionData is 12600(genes) x 181(samples)

expressionData <- gse2553$GSE2553_series_matrix.txt.gz@assayData$exprs 
classes <- gse2553$GSE2553_series_matrix.txt.gz@phenoData@data$description

## create a dataframe [rows: samples, columns: genes(characteristics/features)]
## select only samples having types "Liposarcoma" "Malignant Fibrous Histiocytoma"
## 66 samples are returned : a matrix 66x12601

result <-as.data.frame(cbind(classes, t(expressionData)))
in_data <- result[result$classes %in% c("Diagnosis: Liposarcoma", "Diagnosis: Malignant Fibrous Histiocytoma"), ]

print(paste("in_data - before omit NA vals: ", nrow(in_data), "X", ncol(in_data)))

##
##  omit rows from input data having NA values
##  na.omit removes ROWs on columns having na values
##  so we need to have genes in the ROWs
in_data <- na.omit( as.data.frame(t(in_data)))
print(paste("in_data - after omit NA vals ", nrow(in_data), "X", ncol(in_data)))

in_data <- as.data.frame(t(in_data))
Diagnosis <- as.matrix(in_data$classes)
colnames(Diagnosis) <- c("Diagnosis")
in_data <- in_data[, -1]
print(paste("In Data (after removing classes) ", nrow(in_data), "X", ncol(in_data)))

sample_names <- row.names(in_data)
gene_ids <- as.matrix(colnames(in_data)) 


##  get only numeric
## 
## we need to have in_data as a simple matrix for the function to work

tsne_data <- as.matrix(sapply(in_data, as.numeric))
library(Rtsne)
tsne_out <- Rtsne(tsne_data, perplexity = 20, dim=2) # Run TSNE
plot(tsne_out$Y[, 1] ~ tsne_out$Y[, 2], pch = 20, col = "black")

## apply correlation 
## cor_data matrix is a 8594x8594 matrix
cor_data <- cor(tsne_data)
## save corr data into a tmp triangular matrix  set upper part and diagonal to 0
tmp <- cor_data
tmp[upper.tri(tmp)] <- 0
diag(tmp) <- 0

## exclude columns where correllation > 0.5
selected_data <- in_data[, !apply(tmp, 2, function(x) any(abs(x) > 0.5, na.rm = TRUE))]
selected_data_rows <- nrow(selected_data)
selected_data_cols <- ncol(selected_data)

## selected data is 66x413
print(paste("applying abs(corr) should be <0.5 results in ", selected_data_cols, "genes found"))
print(paste("Selected Data ", nrow(selected_data), "X", ncol(selected_data)))

## apply Rtsne to the selected of the input data (corr <= 0.5)
tsne_out_corr <- Rtsne(selected_data, perplexity = 20, dim=2) # Run TSNE
plot(tsne_out$Y[, 1] ~ tsne_out$Y[, 2], pch = 20, col = "red")




## Thema 4

library(caret)
data <- in_data
data <- lapply(data, as.numeric)
data <- as.data.frame(data)
Diagnosis <- as.factor(Diagnosis)
data <- cbind(in_data, Diagnosis)

inTraining <- caret::createDataPartition(data$Diagnosis, p = .75, list = FALSE)

training <- data[inTraining, ]
testing  <- data[-inTraining, ]

print(levels(training$Diagnosis))

## 10-fold CV repeated 10 times
fitControl <- trainControl( method = "repeatedcv", number = 10, repeats = 10)

names(data)
model_knn1 <- train(training[,1:ncol(training)-1], training$Diagnosis, method='knn', trControl = fitControl)

set.seed(825)

model_knn2 <- train(Diagnosis ~ ., 
                   data = data, 
                   method = "knn", 
                   trControl = fitControl)

#names(getModelInfo())

model_knn <- train(X_train, y_train, method='knn')
prediction_knn<-predict.train(object=model_knn, X_test, type="raw")
table(prediction_knn)
table(prediction_knn, y_test)
confusionMatrix(prediction_knn, as.factor(y_test))
# Train an SVM
model_svm <- train(X_train, y_train, method='svmLinear')
predict_svm<-predict.train(object=model_svm, X_test, type="raw")
table(predict_svm)
table(predict_svm, y_test)
#confusionMatrix(predict_svm,y_test)


## Thema 2
set.seed(1) 
## scale selected data (66x8459)
## 
scaled_data <- sapply(in_data, as.numeric)
scaled_data <- scale(scaled_data)
rownames(scaled_data) <- row.names(in_data)
print(paste("scaled selected data, rows:", nrow(scaled_data), " columns:", ncol(scaled_data)))

print("Find best kmeans cluster separation  k = 1..12 scaled selected data")
scaled_result <- factoextra::fviz_nbclust(scaled_data, kmeans, method = "wss", k.max = 12)
scaled_result
scaled_result$data

## run k-means for scale data & cluster centers
km_selected_scaled <- kmeans(scaled_data, centers = 6, nstart = 25)
## display clusters
clusters_inf_scaled <- factoextra::fviz_cluster(km_selected_scaled, data = scaled_data)
clusters_inf_scaled


##
# Dissimilarity matrix
d <- dist(scaled_data, method = "euclidean")
# Hierarchical clustering using single Linkage
hc1 <- hclust(d, method = "single" )

# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1,  main="Dendrogram Scaled Data - samples 66x8594")

scaled_tsne <- as.data.frame(scale(tsne_out$Y))
rownames(scaled_tsne) <- row.names(in_data)
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

