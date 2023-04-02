########### create random data seq of 5000 b ###################
set.seed(1)
data <-stringi::stri_rand_strings(1, 5000, '[ACGT]')

##### calculate Entropy of the original sequence            ########
##### we could use this instead of Hmax = 2 that assumes    ########
##### equipropable events with prob 0.25 each               ########
H_derived <- 0                                      # init to 0
seq_elements <- c("A", "G", "C", "T")               # we have 4 base noucleotides
distribution_matrix <- matrix(0, 1, 4)              # used to store p(element_occur)*log2(element_occur), elements: G,C,T,A
for (i in 1:length(seq_elements)) {
  distribution_matrix[1, i] = stringr::str_count(toString(data), seq_elements[i])     # get count
  distribution_matrix[1, i] <- distribution_matrix[1, i] / 5000                       # get freq 
  distribution_matrix[1, i] <- distribution_matrix[1, i] * log2(distribution_matrix[1, i])
  H_derived <- H_derived + distribution_matrix[1, i]
}
H_derived <- H_derived *(-1)    ## by definition 
print(H_derived)


no_rows    <-500              # break the initial sequence of 5000 chars
no_columns <- 10              # into 500 substrs each having length 10
data2 <- matrix(nrow = no_rows, ncol = no_columns)
for(i in 1:no_rows) {
  x <- substr(data, (i-1)*10+1, i*10)
  data2[i,]<- unlist(strsplit(x, split = ""))
}
class(toString(data2[, 1]))
################## create pwm ######################################

pwm_matrix <- matrix(0, 4, 10)

rownames(pwm_matrix) <- seq_elements  ### set row names to seq_elements: ("A", "G", "C", "T") we have 4 base noucleotides
colnames(pwm_matrix) <- c("pos1", "pos2", "pos3", "pos4", "pos5", "pos6", "pos7", "pos8", "pos9", "pos10")

for (j in 1:no_columns) {
  for (i in 1:length(seq_elements)) {
    pwm_matrix[i,j] = (stringr::str_count(toString(data2[, j]), seq_elements[i]))/no_rows
  }
}
writeLines("\n")
print(pwm_matrix)
writeLines("\n")
####### export results to excel ##########
library(writexl)
write_xlsx(as.data.frame(pwm_matrix), 'pwm_matrix.xlsx')

####################################### Shannon Entropy ##############################################
## create a 2x10 entropy & info content matrix
## first row: entropy at each pos 1 to 10
## second row: info content at each pos 1 to 10.
######################################################################################################

#### set the table ###
entropy_info_matrix <- matrix(0, nrow = 2, ncol = no_columns)    # entropy and info content matrix
colnames(entropy_info_matrix) <- c("pos1", "pos2", "pos3", "pos4", "pos5", "pos6", "pos7", "pos8", "pos9", "pos10")
rownames(entropy_info_matrix) <- rbind("cell entropy", "info cell content")

## calculate entropy at each pos 1..10 of the motif
for (j in 1:no_columns) {
  for (i in 1:length(seq_elements)) {
    entropy_info_matrix[1, j] = entropy_info_matrix[1, j] + (pwm_matrix[i,j] * log2(pwm_matrix[i,j]))
  }
  entropy_info_matrix[1, j] = (-1) * entropy_info_matrix[1, j]  # by definition 
}


### for info content for pos 1..10 of the motif - we calculate 2 cases ####
###    first:   assuming H initial = Hmax = 2 i.e. euiprobable events with probability 0.25
###   second:   using H_derived above which very close to 2 since we randomly initialized our string 
###
### print also total info content for each case


Hvec <- c(2, H_derived)
for (h in 1:2) {
  dH <- Hvec[h]                # run for each H  
  total_info_content <- 0
  
  for (j in 1:no_columns) {
    entropy_info_matrix[2, j] = dH - entropy_info_matrix[1, j]  # info content
    total_info_content <- total_info_content + entropy_info_matrix[2, j]  
  }   
  writeLines("\n")
  print(paste("Entropy - Info content matrix using H = ", toString(Hvec[h]), sep = ""))
  print(entropy_info_matrix)
  writeLines("\n")
  
  print(paste("Total info content for H = ", toString(Hvec[h]), "is - total = ", toString(total_info_content, sep = "")))
  write_xlsx(as.data.frame(entropy_info_matrix), paste("entropy_info_matrix", toString(h),".xlsx", sep="_"))

}

