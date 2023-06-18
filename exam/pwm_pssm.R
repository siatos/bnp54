print('hopefully copying correctly the [given] sequence ...')
sequence_length <- 50
data <- "AATGCATTGGAATGCTATGCAATGCAATCGAAAGGTCAGGAATGCAATGC"


str(data)
seq_elements <- c("A", "G", "C", "T")   

no_rows    <-10              # break the initial sequence of 50 chars
no_columns <- 5              # into 10 substrs each having length 5
data2 <- matrix(nrow = no_rows, ncol = no_columns)
str(data2)
for(i in 1:no_rows) {
  x <- substr(data, (i-1)*no_columns+1, i*no_columns)
  print(x)
  data2[i,]<- unlist(strsplit(x, split = ""))
}
class(toString(data2[, 1]))
################## create pwm ######################################

pwm_matrix <- matrix(0, 4, no_columns)
print(pwm_matrix)

rownames(pwm_matrix) <- seq_elements  ### set row names to seq_elements: ("A", "G", "C", "T") we have 4 base noucleotides
colnames(pwm_matrix) <- c("pos1", "pos2", "pos3", "pos4", "pos5")

for (j in 1:no_columns) {
  for (i in 1:length(seq_elements)) {
    pwm_matrix[i,j] <- length(grep(seq_elements[i], data2[,j]))/no_rows
  }
}
writeLines("\n")
print(pwm_matrix)
writeLines("\n")

pssm_matrix <- matrix(0, 4, no_columns)
for (j in 1:no_columns) {
  print(j)
  for (i in 1:length(seq_elements)) {
    pssm_matrix[i,j] <- log2(pwm_matrix[i,j]/0.25)
  }
}
print(pssm_matrix)
pssm_matrix[!is.finite(pssm_matrix)] <- 0
print(pssm_matrix)
