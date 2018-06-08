## for all the functions used in SeqQuery




key_extraction = function(mychr, mypos_start, mypos_end, tgf_data)
{
  # mychr = "chr7"
  # mypos_start = 55146606
  # mypos_end = 55146740
  # tgf_data = readRDS("/data/ginny/GENCODE/input_files/tgf_df_CDS.Rds")
  #

  this = tgf_data %>%
    dplyr::filter(chr == mychr) %>%
    dplyr::filter(start_position <= mypos_start) %>%
    dplyr::filter(end_position >= mypos_end)


  return(this)

}







get_nucleotides = function(myseq, mypos_start, mypos_end)
{
  mynucleotides = substr(myseq, mypos_start, mypos_end)
  return(mynucleotides)

}






get_translation = function(excerpts, mypos_start,mypos_end, myseq, dictionary)
{

#   excerpts = this
#   mypos_start = 55146606
#   mypos_end = 55146740
#   myseq = readRDS("/data/ginny/GENCODE/input_files/chromosomes/chr7_seq.Rds")
#   dictionary = readRDS("/data/ginny/GENCODE/input_files/codon_dictionary.Rds")

  # dictionary = triplet_dic


  result =  excerpts %>%
    dplyr::mutate(position_start = mypos_start, position_end = mypos_end)%>%
    dplyr::mutate(query_nucleotide = "",
                  query_codon = "", query_codon_pos = 0, query_aa = "") %>%
    dplyr::select(chr, position_start,position_end, query_nucleotide,
                  query_codon, query_codon_pos, query_aa, everything())


  for(i in 1:nrow(excerpts))
  {

    this_start = excerpts$start_position[i]
    this_end = excerpts$end_position[i]



    segment = get_nucleotides(myseq = myseq,
                              mypos_start = this_start,
                              mypos_end = this_end)

    result$query_nucleotide[i] = substr(segment, (mypos_start-this_start+1), (mypos_end-this_start+1) )

    ls = nchar(segment)

    ef = as.numeric(excerpts$frame[i])

    trans_segment = substr(segment, ef+1, ls)

    ### substitute "T" with "U"

    sub_trans_segment = gsub("T", "U", trans_segment)

    codons = substring(sub_trans_segment,
                       seq(1,nchar(sub_trans_segment),3),
                       seq(3,nchar(sub_trans_segment),3))

    #get_seq = paste0(dictionary[codons], collapse = "")
    #result$translated_seq[i] = get_seq

    ## the below need to quzheng, think over this

    which_codon_start =  floor((mypos_start-this_start-ef)/3+1)   ## I think I need to consider frame here
    which_codon_end = floor((mypos_end-this_start-ef)/3+1)
    result$query_codon[i] = paste(codons[which_codon_start:which_codon_end], collapse = "_")

    start_codon_pos = 0
    end_codon_pos = 0

    if(mypos_start-this_start-ef>0)
    {
      start_codon_pos = (mypos_start-this_start-ef)%%3+1
    }

    if(mypos_end-this_start-ef>0)
    {
      end_codon_pos = (mypos_end-this_start-ef)%%3+1
    }

    result$query_codon_pos[i] = paste0(start_codon_pos, "_",end_codon_pos)

    result$query_aa[i] =  paste(dictionary[codons[which_codon_start:which_codon_end]], collapse = "_")

  }


  return(result)

}



# ### I will build another function to track the changes on the protein 
# 
# get_mutation_result = function(translation_result, ref_nuc, alt_nuc)
# {
#   # if the altered nuc is of the same length as the ref_nuc 
#   # I dont need to extract more nucleotides
#   
#   filter_result = translation_result %>%
#     dplyr::filter(query_nucleotide == ref_nuc)
#   
#   if(nrow(filter_result)>0)
#   {
#     rbindlist(lapply(1:nrow(filter_result), function(x) {
#       alt_nuc_rna = gsub("T","U",filter_result$alt_nuc[x])
#       
#       connect_codon = gsub("_","",filter_result$query_codon[x])
#       
#       two_pos = unlist(strsplit(filter_result$query_codon_pos[x], split = "_", fixed = T))
#       
#       
#      num_codon = 1+lengths(regmatches(filter_result$query_codon[x], gregexpr("_", filter_result$query_codon[x])))
#       
#       sub_start = two_pos[1]
#       sub_end = (num_codon-1)*3+two_pos[2]
#       
#    substr(connect_codon, sub_start, sub_end) <- alt_nuc_rna
#       
#    ### then translate the new codon and you are done.
#    
#    
#       
#    
#    
#     })
#     
#     
#     
#     
#     
#   }
#   
#   
#   
# }
# 
# 
# 




### so only the funciton below will be used by user


#' query function
#'
#' This function allows you to query a chormosome and a position
#' @param query_chr a string indicating the chromosome you want to query
#' @param query_pos_start an integer indicating the starting position you are intersted
#' @param query_pos_end an integer indicating the ending position you are intersted
#' @param proc_tgf_path a string indicating the path to the processed tgf file (provided by the package,"tgf_df_CDS.Rds")
#' @param proc_chr_path a string indicating the path to the directory of chromose fasta files (provided by the package,"chr*_seq.Rds")
#' @param codon_dictionary_path a string indicating the path to the codon directory (provided by the package, "codon_dictionary.Rds")
#' @import dplyr magrittr data.table
#' @export
#' @examples
#' getit = seq_query(query_chr = "chr1",
#'                   query_pos_start = 70000,
#'                   query_pos_end = 70000,
#'                   proc_tgf_path = "/data/ginny/GENCODE",
#'                   proc_chr_path = "/data/ginny/GENCODE/chromosomes",
#'                   codon_dictionary_path = "/data/ginny/GENCODE")



seq_query = function(query_chr, query_pos_start, query_pos_end, proc_tgf_path, proc_chr_path, codon_dictionary_path)
{


  #### I think it is better to have every readin and out only once

  # query_chr = "chr10"
  # query_pos_start = 248336
  # query_pos_end  = 248342
  # proc_tgf_path = "/data/ginny/GENCODE"
  # proc_chr_path = "/data/ginny/GENCODE/chromosomes"
  #

  if(query_pos_start > query_pos_end)
  {
    stop("Invalid position input! ")
  }

  wseq = c(seq(1:22),"X","Y")
  wchr = paste0("chr", wseq)

  if(query_chr %in% wchr == F)
  {
    stop("Invalid chromosome input! ")
  }

    myfilename = paste0(proc_chr_path, "/",query_chr,"_seq.Rds")
    dna_seq = readRDS(myfilename)

  if(query_pos_start<1 | query_pos_end > nchar(dna_seq))
  {
    stop("Invalid position input! ")
  }


  see_nuc = get_nucleotides(myseq = dna_seq,
                            mypos_start = query_pos_start,
                            mypos_end = query_pos_end )

  cat("Query nucleotide(s): ", see_nuc, "\n")




  proc_tgf = readRDS(paste0(proc_tgf_path, "/tgf_df_CDS.Rds"))
  codon_dictionary = readRDS(paste0(codon_dictionary_path,"/codon_dictionary.Rds"))

  see_CDS = key_extraction(mychr = query_chr,
                           mypos_start = query_pos_start,
                           mypos_end = query_pos_end,
                           tgf_data = proc_tgf)

  see_aa = NULL
  if(nrow(see_CDS)>0)
  {
    see_aa = get_translation(excerpts = see_CDS,
                             mypos_start  = query_pos_start,
                             mypos_end = query_pos_end,
                             myseq = dna_seq,
                             dictionary = codon_dictionary)
    cat("Return transcription and translation information!" ,"\n")

  }else{
    cat("Query nucleotide(s) not in CDS!","\n")
  }

  return(see_aa)

}



#
# getit = seq_query(query_chr = "chr1",
#                   query_pos_start = 1419546,
#                   query_pos_end = 1419548,
#                   proc_tgf_path = "/data/ginny/GENCODE",
#                   proc_chr_path = "/data/ginny/GENCODE/chromosomes",
#                   codon_dictionary_path = "/data/ginny/GENCODE")

