####################### create (same) random seq - using same seed always #########################
set.seed(1)
sequence_length <- 5000
data <-stringi::stri_rand_strings(1, sequence_length, '[ACGT]')
class(data)


### First Question #################################################################################
### In the given sequence
###   a. find the distribution of the individual elements
###   b. verify from results if 2nd rule of Chargaff is followed
print("Ans: 1st Q =================================================================================")
####################################################################################################

####################################### Initialize vars #############################################
### For the 4 elements i.e. "A", "G", "T", "C" of the sequence create a 2x4 matrix         #########
### to hold number of appearances i.e. count and frequency                                 #########
####################################################################################################
distribution_matrix <- matrix(0, 2, 4)     # create a 2x4 matrix and initialize to 0
seq_elements <- c("A", "G", "T", "C")
colnames(distribution_matrix) <- seq_elements                            # name the columns
rownames(distribution_matrix) <- rbind("element_count", "element_freq")  # name the rows

## count i.e. occurrence of individual elements i.e. "A", "G", "T", "C" in the sequence data    ###
for (i in 1:length(seq_elements)) {
  distribution_matrix[1, i] = stringr::str_count(toString(data), seq_elements[i])     # get count
  distribution_matrix[2, i] = distribution_matrix[1, i] / sequence_length     # get frequency
}

################################ print results  ####################################################
writeLines("\n")
print(distribution_matrix)
writeLines("\n")
####### export results to excel ##########
library(writexl)
write_xlsx(as.data.frame(distribution_matrix), 'distribution_matrix.xlsx')
####### export results to excel ##########
writeLines("\n")
sprintf("2nd Chargaff law (occ): A~T: % d - % d  G~C: % d - %d", distribution_matrix[1, 1], distribution_matrix[1, 3], distribution_matrix[1,2], distribution_matrix[1, 4])
sprintf("2nd Chargaff law (frq): A~T: % f - % f  G~C: % f - %f", distribution_matrix[2, 1], distribution_matrix[2, 3], distribution_matrix[2,2], distribution_matrix[2, 4])
writeLines("\n")




### Second Question ################################################################################
### split given data sequence into 10 substrings each of 500 chars (500 x 10)
###     find number of occurrences of GC pattern in each of the 10 substrings
print("Ans: 2nd Q  G+C content in 10 substrs of length 500 ========================================")
####################################################################################################

sub_length <- 500
no_of_substrs <-  sequence_length / sub_length
GC_seq_elements <- c("C", "G")

## define a 2x10 matrix to hold
##   first row:  number of "G & C" occurrences per substr of length 500
##   second row: frequency  of "G & C" occurrences per substr of length 500
pattern_distrib_matrix <- matrix(0, 2, no_of_substrs)
colnames(pattern_distrib_matrix) <- c("sseq1", "sseq2", "sseq3", "sseq4", "sseq5", "sseq", "sseq7", "sseq8", "sseq9", "sseq10")
rownames(pattern_distrib_matrix) <- rbind("number of G+C content", "relative freq of G+C content")

for (i in 1:no_of_substrs) {                 # instead of a single loop with "i in seq format"
   start_str <- (i-1)*sub_length+1           # i.e. 1, 501, ..., 4501  etc - calc start of substr
   end_str   <- start_str + sub_length -1    # i.e  500, 1000, ..., 4500 etc - calc end of substr
   for (j in 1:length(GC_seq_elements)) {
      pattern_distrib_matrix[1, i] <- pattern_distrib_matrix[1, i] + stringr::str_count(substr(data, start_str, end_str), GC_seq_elements[j])
   }
}
for (i in 1:no_of_substrs){
  pattern_distrib_matrix[2, i] <- pattern_distrib_matrix[1, i] / sub_length     # get frequency
}

writeLines("\n")
print(pattern_distrib_matrix)
writeLines("\n")
####### export results to excel ##########
write_xlsx(as.data.frame(pattern_distrib_matrix), 'pattern_distrib_matrix.xlsx')
####### export results to excel ##########

###   (kind of of) simple scatter plot
xtick<-seq(1, 10, by=1)
plot(pattern_distrib_matrix[2,], col = "blue", pch = 19, main = "visualize G+C content", ylab = "G+C freq", xlab = "id pos in subseq")
axis(1, xtick)
meanv <- mean(pattern_distrib_matrix[2,])
stdv  <- sd(pattern_distrib_matrix[2,])
abline(h=meanv, lwd=1.5, col="red")

abline(h=meanv-stdv, lwd=1.5, col="green")
abline(h=meanv+stdv, lwd=1.5, col="green")
