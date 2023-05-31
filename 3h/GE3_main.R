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




# Thema 3
# need to install
#BiocManager::install("randomForest")
# Set RFE control

data <- in_data
data[] <- lapply(data, as.numeric)  ## convert to numeric
data <- as.data.frame(data)
Diagnosis <- as.factor(Diagnosis)
data <- cbind(data, Diagnosis)

ctrl = rfeControl(functions = rfFuncs, # "rfFuncs" are built-in to caret
                  method = "repeatedcv", repeats = 10,
                  saveDetails = TRUE)
# By using rfFuncs, caret will use a random forest to evaluate the usefulness of a feature.
# Set a sequence sizes to search
sizes=c(50, 75, 100, 150, 200, 300, 400, 500, 600, 1000)

# Use caret's rfe function to fit RF models to these different feature spaces
rfeResults = rfe(x = data[,1:ncol(data)-1], y = data[,ncol(data)],
                 sizes = sizes,
                 rfeControl = ctrl)


#rfeResults = rfe(x = select(nonColinearData, -Depressiv), y = nonColinearData$Depressiv,
#                 sizes = sizes,
#                 rfeControl = ctrl)


rfeResults$results
ggplot(data = rfeResults, metric = "Accuracy") + theme_bw()
ggplot(data = rfeResults, metric = "Kappa") + theme_bw()

X <- rfeResults$variables
Y <- X[X$Variables == 1000, ]
Y <- Y[Y$Resample == "Fold01.Rep01", ]
final_data <- in_data[, Y$var]



#f = rfeResults.get_support(1) #the most important features
#X = data[data.columns[f]]


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
fitControl <- trainControl( method = "repeatedcv", number = 10, repeats = 10)

set.seed(1)
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
        

        

####################################################


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



