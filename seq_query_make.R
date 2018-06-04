library("devtools")
devtools::install_github("klutometis/roxygen")
library(roxygen2)

setwd("/Users/ginny")
create("SeqQuery")


setwd("/Users/ginny/SeqQuery")
document()


setwd("/Users/ginny")
install("SeqQuery")




library("devtools")
devtools::install_github("ginnyintifa/SeqQuery")
library(SeqQuery)



# 
# 
# getit = seq_query(query_chr = "chr1",
#                   query_pos_start = 1419546,
#                   query_pos_end = 1419548,
#                   proc_tgf_path = "/data/ginny/GENCODE/input_files",
#                   proc_chr_path = "/data/ginny/GENCODE/input_files/chromosomes",
#                   codon_dictionary_path = "/data/ginny/GENCODE/input_files")
# 


