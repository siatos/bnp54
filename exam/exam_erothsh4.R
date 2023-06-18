## full path for the needed input file
## needs to be (can be) adjusted accordingly
## here we use "root" folder c:/

exercise_txt_path <- "C:/data_bnp54.txt"
library(readr)
## read files input and class data in a data frame
in_data <- read.table(exercise_txt_path, sep="\t", header=TRUE)
classes_txt_path <- "C:/class_bnp54.txt"
in_classes <- read.table(exercise_txt_path, sep="\t", header=TRUE)
str(in_data)
colnames <- colnames(in_data)
str(in_classes)

# Load libraries
library(cluster)
library(factoextra)
library(rpart)

in_data1 <- in_data[, 1:10]
tree_model1 <- rpart( ~ ., data = in_data1, method = "in_classes")


in_data2 <- in_data[, 91:100]
tree_model1 <- rpart( ~ ., data = in_data1, method = "in_classes")


library(caret)
inTraining <- caret::createDataPartition(data$Diagnosis, p = .75, list = FALSE)
training <- data[inTraining, ]
testing  <- data[-inTraining, ]
