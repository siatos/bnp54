########################################################################################
##
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

## save Diagnosis as column
Diagnosis <- as.matrix(in_data$classes)
colnames(Diagnosis) <- c("Diagnosis")
Diagnosis <- as.factor(Diagnosis)
in_data <- in_data[, -1]
print(paste("In Data (after removing classes) ", nrow(in_data), "X", ncol(in_data)))

## save gene_ids
sample_names <- row.names(in_data)
gene_ids <- as.matrix(colnames(in_data)) 

##  get only numeric
## 
## we need to have in_data as a simple matrix for the tsne function to work
set.seed(4300)
tsne_data <- as.matrix(sapply(in_data, as.numeric))
library(Rtsne)
tsne_out <- Rtsne(tsne_data, perplexity = 20, dim=2) # Run TSNE
tsne_plot <- data.frame(x = tsne_out$Y[ ,1], y = tsne_out$Y[, 2], col = Diagnosis) 
library(ggplot2)
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col)) + theme(legend.position="right") + scale_color_manual(values = c("blue", "red"))


## apply correlation 
## cor_data matrix is a 8594x8594 matrix
cor_data <- cor(tsne_data)

## save corr data into a tmp triangular matrix  (matrix is symmetrical over the diagaonal) 
## set upper part and diagonal to 0
tmp <- cor_data
tmp[upper.tri(tmp)] <- 0
diag(tmp) <- 0

## exclude columns where abs(correllation) > 0.5
selected_data <- in_data[, !apply(tmp, 2, function(x) any(abs(x) > 0.5, na.rm = TRUE))]
selected_data_rows <- nrow(selected_data)
selected_data_cols <- ncol(selected_data)

## selected data is 66x413
print(paste("applying abs(corr) should be <0.5 results in ", selected_data_cols, "genes found"))
print(paste("Selected Data ", nrow(selected_data), "X", ncol(selected_data)))
set.seed(4300)
## apply Rtsne to the selected of the input data (corr <= 0.5)
tsne_out_corr <- Rtsne(selected_data, perplexity = 20, dim=2) # Run TSNE
tsne_plot <- data.frame(x = tsne_out_corr$Y[,1], y = tsne_out_corr$Y[,2], col = Diagnosis) 
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col)) + theme(legend.position="right") + scale_color_manual(values = c("blue", "red"))


########################################################################################


########################################################################################
##
## Thema 2
set.seed(1) 
## scale selected data (66x8594)
scaled_data <- sapply(in_data, as.numeric)
scaled_data <- scale(scaled_data)
rownames(scaled_data) <- row.names(in_data)
print(paste("scaled selected data, rows:", nrow(scaled_data), " columns:", ncol(scaled_data)))

##
## apply k-means 
print("Find best kmeans cluster separation  k = 1..12 on scaled selected data")
scaled_result <- factoextra::fviz_nbclust(scaled_data, kmeans, method = "wss", k.max = 12)
scaled_result
scaled_result$data

## run kmeans for scaled data & cluster centers
## select k=7 as the optimal solution - previous elbow plot 
##
km_res <- kmeans(scaled_data, centers = 7, nstart = 25)
## display clusters
kmclusters <- factoextra::fviz_cluster(km_res, data = scaled_data)
kmclusters

##
## apply agglomerative hcl
# 
library("cluster")
set.seed(1)

x<- sapply(in_data, as.numeric) #x <- scaled_data
rownames(x) <-row.names(in_data)
res.agnes <- agnes(x, diss = FALSE,stand = TRUE, metric = "euclidean",  method = "single")   
x.grps <- cutree(res.agnes, 2:12)
res_m <- matrix(0, nrow = 11, ncol = 2)
for (i in 1:11) {
  x.SS <- aggregate(x, by=list(x.grps[, i]), function(x) sum(scale(x, scale=FALSE)^2))
  SS <- rowSums(x.SS[, -1])   # Sum of squares for each cluster
  TotalSS <- sum(x.SS[, -1])  # wss 
  res_m[i, 1] <- i+1
  res_m[i, 2] <- TotalSS
}
res_m
library(ggplot2)
plot_data <- data.frame(X=as.factor(res_m[,1]), Y=res_m[,2])
ggplot(plot_data, aes(x = X, y = Y, group = 1)) +
  geom_line() +
  geom_point()+
  xlab("no of Clusters") +
  ylab("wss") 
hc_res <- factoextra::hcut(x, k = 7, hc_func = "agnes", hc_method = "single", hc_metric = "euclidean")
kcolors <- c("blue", "yellow", "magenta", "green", "black", "red", "yellow")
#kcolors <- c("blue", "yellow", "magenta", "gray", "green", "black", "red", "pink", "orange", "purple")

##
## display dendrogram
factoextra::fviz_dend(hc_res, k = 7, # Cut in 7 groups
                      cex = 0.5, # label size
                      k_colors = kcolors,
                      color_labels_by_k = TRUE, # color labels by groups
                      rect = TRUE, # Add rectangle around groups
                      rect_fill = TRUE)


p1 <- factoextra::fviz_nbclust(scaled_data, FUN = factoextra::hcut, method = "wss", k.max = 11) + ggtitle("(A) Elbow method")
p1

## apply k-means on tsne results
##
rownames(tsne_out$Y) <- row.names(in_data)
tsne_result <- factoextra::fviz_nbclust(tsne_out$Y, kmeans, method = "wss", k.max = 12)
tsne_result
km_tsne <- kmeans(tsne_out$Y, centers = 4, nstart = 25)
## display clusters
clusters_tsne <- factoextra::fviz_cluster(km_tsne, data = as.data.frame(tsne_out$Y))
clusters_tsne

##
########################################################################################


########################################################################################
##
## Thema 3
## Set RFE control
set.seed(4300) #1200
library(caret)
data <- in_data
data[] <- lapply(data, as.numeric)  ## convert to numeric
data <- as.data.frame(data)
# Diagnosis <- as.factor(Diagnosis)
data <- cbind(data, Diagnosis)

ctrl = rfeControl(functions = rfFuncs, # "rfFuncs" are built-in to caret
                  method = "repeatedcv", repeats = 5,
                  saveDetails = TRUE)

# By using rfFuncs, caret will use a random forest to evaluate features
# Set a sequence sizes to search

#sizes=c(50, 75, 100, 150, 200, 300, 400, 500, 600, 1000)
sizes=c(50, 200, 500, 1000)

# Use caret's rfe function to fit RF models to these different feature spaces
rfeResults = rfe(x = data[,1:ncol(data)-1], y = data[,ncol(data)],
                 sizes = sizes,
                 rfeControl = ctrl)

rfeResults$results
#dev.off()
ggplot(data = rfeResults, metric = "Accuracy") + theme_bw()
ggplot(data = rfeResults, metric = "Kappa") + theme_bw()

X <- rfeResults$variables
Y <- X[X$Variables == 1000, ]
Y_1000 <- aggregate(Y[, c("Overall")], list(Y$var), mean)
Y_1000 <- order(Y_1000[, c("x")], decreasing = TRUE)[1:1000]
selected_1000 <- in_data[, Y_1000]

tsne_1000 <- Rtsne(selected_1000, perplexity = 20, dim=2) # Run TSNE
tsne_plot <- data.frame(x = tsne_1000$Y[,1], y = tsne_1000$Y[,2], col = Diagnosis) 
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col)) + theme(legend.position="right") + scale_color_manual(values = c("blue", "red"))


########################################################################################



########################################################################################
##
## Thema 4

library(caret)
data <- in_data
data[] <- lapply(data, as.numeric)  ## convert to numeric
data <- as.data.frame(data)
Diagnosis <- as.factor(Diagnosis)
data <- cbind(data, Diagnosis)

inTraining <- caret::createDataPartition(data$Diagnosis, p = .75, list = FALSE)

training <- data[inTraining, ]
testing  <- data[-inTraining, ]

print(levels(training$Diagnosis))

x_train <- training[,1:ncol(training)-1]  # data[idxs,1:4]
y_train <- training[ , ncol(training)]    # training$Diagmosis  
x_test  <- testing[, 1:ncol(testing)-1]    # data[-idxs,1:4]
y_test  <- testing[, ncol(testing)]       # testing$Diagmosis   

# prepare the matrix for results
algorithm.knn <- c(0, 0, 0)
algorithm.svm <- c(0, 0, 0)
algorithm.rf  <- c(0, 0, 0)

## 10-fold CV repeated 10 times
## fitControl <- trainControl( method = "repeatedcv", number = 10, repeats = 10)
## 5-fold CV repeated 5 times
fitControl <- trainControl( method = "repeatedcv", number = 5, repeats = 5)


set.seed(121)
model_knn <- train(x_train, y_train, method='knn', trControl = fitControl)
prediction_knn <- predict.train(object=model_knn, x_test, type="raw")
table(prediction_knn, y_test)
knn_result  <- confusionMatrix(prediction_knn, as.factor(y_test))
knn_result$byClass
algorithm.knn <- c(knn_result[["overall"]][["Accuracy"]], knn_result$byClass['Precision'], knn_result$byClass['Recall'])

model_svm <- train(x_train, y_train, method='svmLinear', trControl = fitControl)
prediction_svm <- predict.train(object=model_svm, x_test, type="raw")
table(prediction_svm, y_test)
svm_result  <- confusionMatrix(prediction_svm, as.factor(y_test))
svm_result$byClass
algorithm.svm <- c(svm_result[["overall"]][["Accuracy"]], svm_result$byClass['Precision'], svm_result$byClass['Recall'])

model_rf <- train(x_train, y_train, method='rf', trControl = fitControl)
prediction_rf <- predict.train(object=model_rf, x_test, type="raw")
table(prediction_rf, y_test)
rf_result  <- confusionMatrix(prediction_rf, as.factor(y_test))
rf_result$byClass
algorithm.rf <- c(rf_result[["overall"]][["Accuracy"]], rf_result$byClass['Precision'], rf_result$byClass['Recall'])

algorithms <- rbind(algorithm.knn, algorithm.svm, algorithm.rf)
#algorithms <- as.data.frame(algorithms)
colnames(algorithms) <- c("", "", "")
algorithms

boxplot(algorithms, 
        beside=TRUE, 
        col = c("yellow", "blue", "red"), 
        main = "Algorithm metrics", 
        xlab="Algorithms", 
        ylab="Metrics Values")
        legend("top", legend = c("Accuracy", "Precision", "Recall"), fill = c("yellow", "blue", "red"))
        

########################################################################################





