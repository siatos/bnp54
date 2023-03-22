####################### create (same) random seq - using same seed always #########################
set.seed(1)
sequence_length <- 5000
data <-stringi::stri_rand_strings(1, sequence_length, '[ACGT]')
class(data)

#################################### load necessary libs ###########################################
library(NLP)         # only in order to avoid explicit reference using the package name
library(stringr)
library(stringi)
####################################################################################################

# calc_occurence <- function(strsequence, strpattern)
# {
#     len <- stri_length(strsequence)
#     occ <- str_count(strsequence, strpattern)
#     return(list(occ, occ/len))
# }

### First Question #################################################################################
### In the given sequence
###   a. find the distribution of the individual elements
###   b. verify from results if 2nd rule of Chargaff is followed
print("Ans: 1st Q =================================================================================")
####################################################################################################

####################################### Initialze vars #############################################
### For the 4 elements i.e. "A", "G", "T", "C" of the sequence create a 2x4 matrix         #########
### to hold number of appearances i.e. count and frequency                                 #########
####################################################################################################
distribution_matrix <- matrix(0, 2, 4)     # create a 2x4 matrix and initialize to 0
seq_elements <- c("A", "G", "T", "C")
colnames(distribution_matrix) <- seq_elements                            # name the columns
rownames(distribution_matrix) <- rbind("element_count", "element_freq")  # name the rows

## count i.e. occurence of individual elmements i.e. "A", "G", "T", "C" in the sequence data    ###
for (i in 1:length(seq_elements)) {
  distribution_matrix[1, i] = str_count(as.String(data), seq_elements[i])     # get count
  distribution_matrix[2, i] = distribution_matrix[1, i] / sequence_length     # get frequency
}

################################ print results  ####################################################
writeLines("\n")
print(distribution_matrix)
writeLines("\n")
sprintf("2nd Chargaff law (occ): A~T: % d - % d  G~C: % d - %d", distribution_matrix[1, 1], distribution_matrix[1, 3], distribution_matrix[1,2], distribution_matrix[1, 4])
sprintf("2nd Chargaff law (frq): A~T: % f - % f  G~C: % f - %f", distribution_matrix[2, 1], distribution_matrix[2, 3], distribution_matrix[2,2], distribution_matrix[2, 4])
writeLines("\n")

### Second Question ################################################################################
### split given data sequence into 10 substrings each of 500 chars (500 x 10)
###     find number of occurence of GC pattern in each of the 10 substrings
print("Ans: 2nd Q  G+C content in 10 substrs of length 500 ========================================")
####################################################################################################

sub_length <- 500
no_of_substrs <-  sequence_length / sub_length
GC_seq_elements <- c("C", "G")

## define a 2x10 matrix to hold
##   first row:  number of "G & "C" occurences per substr of length 500
##   second row: frequency  of "G & "C" occurences per substr of length 500
pattern_distrib_matrix <- matrix(0, 2, no_of_substrs)

for (i in 1:no_of_substrs) {                 # instead of a single loop with "i in seq format"
   start_str <- (i-1)*sub_length+1           # i.e. 1, 501, ..., 4501  etc - calc start of substr
   end_str   <- start_str + sub_length -1    # i.e  500, 1000, ..., 4500 etc - calc end of substr
   for (j in 1:length(GC_seq_elements)) {
      pattern_distrib_matrix[1, i] <- pattern_distrib_matrix[1, i] + str_count(substr(data, start_str, end_str), GC_seq_elements[j])
   }
}
for (i in 1:no_of_substrs){
  pattern_distrib_matrix[2, i] <- pattern_distrib_matrix[1, i] / sub_length     # get frequency
}

writeLines("\n")
pattern_distrib_matrix

###   (kind of of) simple scatter plot
x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
plot(pattern_distrib_matrix[2,])