Gene1 <- c(1, 2, 3, 4, 5, 6)
Gene2 <- c(0, 2, 4, 6, 8, 10)
Gene3 <- c(0,1,2,0,1,2)
genes_matrix <- rbind(Gene1, Gene2, Gene3)
rownames(genes_matrix) <- c("Gene1", "Gene2", "Gene3")
colnames(genes_matrix) <- c("Cell1", "Cell2", "Cell3", "Cell4", "Cell5", "Cell6")

data2 <- matrix(nrow = 3, ncol = 11)
for (i in 1:3) {
  Gene_mean <- mean(genes_matrix[i,])
  Gene_max <- max(genes_matrix[i,])
  Gene_min <- min(genes_matrix[i,])
  Gene_quantile_25 <-quantile(genes_matrix[i,], probs = c(0.25))
  Gene_quantile_75 <-quantile(genes_matrix[i,], probs = c(0.75))
  data2[i,] <- c(genes_matrix[i,], Gene_mean, Gene_max, Gene_min, Gene_quantile_25, Gene_quantile_75)
}

rownames(data2) <- c("Gene1", "Gene2", "Gene3")
colnames(data2) <- c("Cell1", "Cell2", "Cell3", "Cell4", "Cell5", "Cell6", "Gene_mean", "Gene_max", "Gen_min", "Gene_quantile_25", "Gene_quantile_75")

row.names(data2)

##### Print the data  #####
print("================================================================")
print("####################### Genes table ############################")
print(genes_matrix)
print("####################### Updated Genes Table #####################")
print("Updated Genes Table ")
print(data2)


##### Create the boxplot  #####
boxplot.matrix(t(genes_matrix), xlab="Genes ", ylab="values")

##### Calculate prop table i.e. contribution of each Gene for every Cell 1..6  #####
print("######################### Prop table ############################")
prop.table(genes_matrix, margin = 2)

##### Creting the barplot : setting colors labels & legend #####
barplot(prop.table(genes_matrix, margin = 2), beside=TRUE, col = c("yellow", "blue", "red"), main = "Gene Fequencies in each Cell", ylim=c(0.0,1.0), xlab="Gene participation on each Cell", ylab="Gene contribution")
legend("top", legend = c("Gene1", "Gene2", "Gene3"), fill = c("yellow", "blue", "red"))

##### Creating the hist :  #####
v <- prop.table(genes_matrix, margin = 2)
v<-prop.table(genes_matrix, margin=2)
hist(v[,1], breaks=50, col=rgb(1, 1,0, 0.4), ylab = "Number of genes", xlab="Total Cell Expression", main="Cell Histogram")  # Cell1  yellow
hist(v[,2], breaks=20, col=rgb(0, 0,1, 0.2), add=TRUE)           # Cell2  blue
hist(v[,3], breaks=20, col=rgb(1, 0,0, 0.2), add=TRUE)           # Cell3  red
hist(v[,4], breaks=20, col=rgb(1, 1,1, 0.2), add=TRUE)           # Cell4  white
hist(v[,5], breaks=20, col=rgb(1, 0,1, 0.2), add=TRUE)           # Cell5  magenta
hist(v[,6], breaks=50, col=rgb(0, 1,0, 0.2), add=TRUE)           # Cell6  green

legend('top', c('Cell1', 'Cell2', 'Cell3', 'Cell4', 'Cell5', 'Cell6'), fill=c('yellow', 'blue', 'red', 'white', 'magenta', 'green'))