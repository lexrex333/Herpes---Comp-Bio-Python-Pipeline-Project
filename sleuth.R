#Sleuth wrapper - Python Pipeline -- taken from slides

#FIRST setting up the directory that all the output files should go to
output_dir = "/home/2025/aavalos4/PythonPipeline_Lexi_Avalos" #change depending on your path 

#move into the directory from now on
setwd(output_dir)

#Load the package
library(sleuth)

#reading in the table I have of samples and kallisto output
directory = file.path(getwd(), "sleuth_input.txt")
stab = read.table(directory, header = TRUE)

#initializing sleuth object 
so = sleuth_prep(stab)

#fitting a model comparing the conditions
so = sleuth_fit(so, ~condition, 'full')

#fitting the reduced model to compare in the likelihood ratio test
so = sleuth_fit(so, ~1, 'reduced')

#performing the likelihood ratio test for differential expression between conditions
so = sleuth_lrt(so, 'reduced', 'full')

#loading the dplyr package for data.frame filtering
library(dplyr)

#extracting the test results from the sleuth object
sleuth_table = sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

#filtering by the most important results (FDR/qval <0.05) and sort by pval
sleuth_significant = dplyr::filter(sleuth_table, qval < 0.05) %>% dplyr::arrange(pval)


#checking - printing top 10 transcripts
#head(dplyr::select(sleuth_significant, target_id, pval, qval), n =10)

#saving and writing significant results
write.table(sleuth_significant[, c("target_id", "test_stat", "pval", "qval")],
                    file= file.path(getwd(), "sleuth_output.tsv"), 
                    sep = "\t", quote = FALSE, row.names = FALSE)