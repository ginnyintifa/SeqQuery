

## 1. Installation

`SeqQuery` can be downloaded and installed in R as follows. As a prerequisite, `devtools` must be installed:

```{r, eval = F}
install.packages("devtools")
```

Next, install `SeqQuery`:

```{r, eval = F}

library("devtools")
devtools::install_github("ginnyintifa/SeqQuery")
library(SeqQuery)

```
## 2. Function

### Query a location/region on a chromosome.

Function ```seq_query``` takes the chromosome and the location range of interest as input parameters. It will return the nucleotide sequence as the first output. If the region of interest happens to be with in CDS (coding DNA sequence), the function will return the amino acids coded by the nucleotides and some other sequence information.


## 3. Input files

### SeqQuery provided feature files

Several files need to be downloaded from the following [WEBSITE](http://137.132.97.109:59739/CSSB_LAB/) before running `SeqQuery`:


* 1 mRNA codon book in **codon_dictionary.Rds**. 
 
* 2 Processed GTF file for information on all the CDS in **tgf_df_CDS.Rds**.
 
* 3 DNA sequences of 24 chromosomes in the folder **chormosomes**, each file **chr*_seq.Rds** contains sequence for each chormosome.
 
The original GTF file and sequence fasta file are downloaded from GENCODE: https://www.gencodegenes.org/releases/current.html
 


## 4. Output

Output from ```seq_query``` is simple. First the nucleotide sequence will be printed. If the queried region happens to be in CDS, then the amino acids coded by the nucleotide sequence along with other transcription and translation information will be the return value of the function.


## 5. Example script

```
getit = seq_query(query_chr = "chr1",
                  query_pos_start = 1419546,
                  query_pos_end = 1419548,
                  proc_tgf_path = "/data/ginny/GENCODE",
                  proc_chr_path = "/data/ginny/GENCODE/chromosomes",
                  codon_dictionary_path = "/data/ginny/GENCODE")
                  
getit                             
```


Output of the above calling:

```
Query nucleotide:  TGT 
Return transcription and translation information! 
> getit
   chr position_start position_end query_nucleotide query_codon query_codon_pos
1 chr1        1419546      1419548              TGT         UGU             2_1
2 chr1        1419546      1419548              TGT         UGU             1_3
3 chr1        1419546      1419548              TGT         AAU             3_2
4 chr1        1419546      1419548              TGT         AAU             3_2
5 chr1        1419546      1419548              TGT         UGU             1_3
  query_aa  source type start_position end_position score strand frame
1        C  HAVANA  CDS        1419545      1419549     .      -     2
2        C  HAVANA  CDS        1419237      1419549     .      -     1
3        N ENSEMBL  CDS        1419103      1419549     .      -     0
4        N  HAVANA  CDS        1419103      1419549     .      -     0
5        C  HAVANA  CDS        1419237      1419549     .      -     1
            gene_id     transcript_id      gene_type gene_name transcript_type
1 ENSG00000235098.8 ENST00000520296.5 protein_coding   ANKRD65  protein_coding
2 ENSG00000235098.8 ENST00000427211.3 protein_coding   ANKRD65  protein_coding
3 ENSG00000235098.8 ENST00000537107.5 protein_coding   ANKRD65  protein_coding
4 ENSG00000235098.8 ENST00000454272.2 protein_coding   ANKRD65  protein_coding
5 ENSG00000235098.8 ENST00000442470.1 protein_coding   ANKRD65  protein_coding
  transcript_name exon_number           exon_id level        protein_id
1     ANKRD65-204           3 ENSE00003617803.1     1 ENSP00000429035.1
2     ANKRD65-201           3 ENSE00003635074.1     2 ENSP00000428419.1
3     ANKRD65-205           4 ENSE00002229397.2     3 ENSP00000445688.1
4     ANKRD65-203           3 ENSE00001784050.1     2 ENSP00000482314.1
5     ANKRD65-202           2 ENSE00001641431.1     2 ENSP00000428201.1
  havana_gene    havana_transcript
1        CCDS          CCDS57963.1
2 CCDS57962.1 OTTHUMG00000002911.4
3        CCDS          CCDS55558.1
4        CCDS          CCDS55558.1
5        CCDS          CCDS57962.1

```









