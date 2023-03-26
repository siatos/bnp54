########### create random data seq of 5000 b ###################
set.seed(1)
data <-stringi::stri_rand_strings(1, 5000, '[ACGT]')

no_rows    <-500              # break the initial sequence of 5000 chars
no_columns <- 10              # into 500 substrs each having length 10
data2 <- matrix(nrow = no_rows, ncol = no_columns)
for(i in 1:no_rows) {
  x <- substr(data, (i-1)*10+1, i*10)
  data2[i,]<- unlist(strsplit(x, split = ""))
}
class(NLP::as.String(data2[, 1]))
################## create pwm ######################################

pwm_matrix <- matrix(0, 4, 10)
seq_elements <- c("A", "G", "C", "T")               # we have 4 base noucleotides
rownames(pwm_matrix) <- seq_elements

for (j in 1:no_columns) {
  for (i in 1:length(seq_elements)) {
    pwm_matrix[i,j] = (stringr::str_count(NLP::as.String(data2[, j]), seq_elements[i]))/no_rows
  }
}
writeLines("\n")
print(pwm_matrix)
writeLines("\n")

####################################### Shannon Entropy ##############################################
## create a 2x10 entropy & info content matrix
## first row: entropy at each pos 1 to 10
## second row: info content at each pos 1 to 0.
######################################################################################################
entropy_info_matrix <- matrix(0, nrow = 2, ncol = no_columns)    # entropy and info content matrix
colnames(entropy_info_matrix) <- c("pos1", "pos2", "pos3", "pos4", "pos5", "pos6", "pos7", "pos8", "pos9", "pos10")
maxH <- 2    # for 4 noucleotides assuming equiproble events
for (j in 1:no_columns) {
  # entropy_info_matrix[1, j] <- 0
  for (i in 1:length(seq_elements)) {
    entropy_info_matrix[1, j] = entropy_info_matrix[1, j] + (pwm_matrix[i,j] * log2(pwm_matrix[i,j]))
  }
  entropy_info_matrix[1, j] = (-1) * entropy_info_matrix[1, j]  # by definition - we need to multiply by (-1) to get > 0
  entropy_info_matrix[2, j] = maxH - entropy_info_matrix[1, j]  # almost equiprobable events have info content close to 0
}
writeLines("\n")
print("Entropy - Info content matrix")
print(entropy_info_matrix)
writeLines("\n")
