# in silico translation 

## data processing script
library(dplyr)
library(magrittr)
library(data.table)

tgf = readLines("/data/ginny/GENCODE/gencode.v27.basic.annotation.gtf")

### only extract the lines with "protein_coding"


tgf_pc = grep("protein_coding", tgf, value = T)


### make the file tab delimited


tgf_df = rbindlist(lapply(1:length(tgf_pc), function(x) {
  
  if(x%%100000 == 0)
    cat(x,"\n")
  
  ss = strsplit(tgf_pc[x], split = "\t", fixed = T)
  
  dss = data.frame(t(unlist(ss)), stringsAsFactors = F)
  
  return(dss)
  
}))




colnames(tgf_df) = c("chr","source","type", "start_position","end_position",
                     "score","strand", "frame", "description")



### only get CDS

tgf_df_CDS = tgf_df %>%
  dplyr::filter(type == "CDS")


parsed_desc = rbindlist(lapply(1:nrow(tgf_df_CDS),function(x) {
  
  this_desc = tgf_df_CDS$description[x]
  
  getit = unlist(strsplit(this_desc, split = ";", fixed = T))
  
  gene_id  = unlist(strsplit(getit[1], split = "\"",fixed = T))[2]
  transcript_id = unlist(strsplit(getit[2], split = "\"",fixed = T))[2]
  gene_type  = unlist(strsplit(getit[3], split = "\"",fixed = T))[2]
  gene_name = unlist(strsplit(getit[4], split = "\"",fixed = T))[2]
  transcript_type = unlist(strsplit(getit[5], split = "\"",fixed = T))[2]
  transcript_name = unlist(strsplit(getit[6], split = "\"",fixed = T))[2]
  exon_number = as.integer(unlist(strsplit(getit[7], split = " ",fixed = T))[3])
  exon_id  = unlist(strsplit(getit[8], split = "\"",fixed = T))[2]
  level = as.integer(unlist(strsplit(getit[9], split = " ",fixed = T))[3])
  protein_id = unlist(strsplit(getit[10], split = "\"",fixed = T))[2]
  havana_gene = unlist(strsplit(getit[14], split = "\"",fixed = T))[2]
  havana_transcript = unlist(strsplit(getit[15], split = "\"",fixed = T))[2]
  
  get_df = data.frame(gene_id, transcript_id, gene_type, gene_name,
                      transcript_type, transcript_name, exon_number, exon_id,
                      level, protein_id, havana_gene, havana_transcript, stringsAsFactors = F)
  
  if(x %% 10000 == 0)
    cat(x, "\n")
  
  return(get_df)
} ))



tgf_df_CDS$start_position = as.integer(tgf_df_CDS$start_position)
tgf_df_CDS$end_position = as.integer(tgf_df_CDS$end_position)


tgf_df_CDS_desc = cbind(tgf_df_CDS, parsed_desc)

tgf_df_CDS_desc = tgf_df_CDS_desc[,-9]

saveRDS(tgf_df_CDS_desc,file = "tgf_df_CDS.Rds")


aa = c("F","F","L","L",
       "S","S","S","S",
       "Y","Y","STOP","STOP",
       "C","C","STOP","W",
       "L","L","L","L",
       "P","P","P","P",
       "H","H","Q","Q",
       "R","R","R","R",
       "I","I","I","M",
       "T","T","T","T",
       "N","N","K","K",
       "S","S","R","R",
       "V","V","V","V",
       "A","A","A","A",
       "D","D","E","E",
       "G","G","G","G")


triple = c("UUU","UUC", "UUA", "UUG",
           "UCU","UCC","UCA","UCG",
           "UAU","UAC","UAA","UAG",
           "UGU","UGC","UGA","UGG",
           "CUU","CUC","CUA","CUG",
           "CCU","CCC","CCA","CCG",
           "CAU","CAC","CAA","CAG",
           "CGU","CGC","CGA","CGG",
           "AUU","AUC","AUA","AUG",
           "ACU","ACC","ACA","ACG",
           "AAU","AAC","AAA","AAG",
           "AGU","AGC","AGA","AGG",
           "GUU","GUC","GUA","GUG",
           "GCU","GCC","GCA","GCG",
           "GAU","GAC","GAA","GAG",
           "GGU","GGC","GGA","GGG")

triplet_dic = structure(aa, names = triple)


saveRDS(triplet_dic, file = "codon_dictionary.Rds")




wf = readLines("all_file/GRCh38.p10.genome.fa")

sig = grepl(">", wf)
starts = which(sig==T)


for(i in 1:24)
{
  startline = starts[i]+1
  endline = starts[i+1]-1
  seq = paste0(wf[startline:endline], collapse = "")
  cat(i, "\n")
  
  saveRDS(seq, file = paste0("chr", i,"_seq.Rds"))
  
}

