require(stringr)
require(plyr)
require(dplyr)
require(data.table)
#
# [0.1] take care about the input file names
#
wfolder <- "/work2/project/sigenae/Project_RNA_Refinery/pipelineOut_embryo/"
contigs_llist <- paste(wfolder, "contigs_after_refseqPEP.llist", sep = "")
blastn_DNA <- paste(wfolder, "contigs_after_CDS_CDNA.PEP.DNA_blasted.table", sep = "")
#
# [0.2] read the input
#
llist <- read.table(
    contigs_llist, header = FALSE, sep = " ", quote = "",
    stringsAsFactors = FALSE,comment.char = "",
    colClasses = c("character","integer"))
# > dim(llist)
# [1] 473920      2
# 
blast_DNA <- fread(blastn_DNA)
# > dim(blast_DNA)
# [1] 264400      9
#
setnames(blast_DNA, c("query","hit","desc","identity","bits","qlen","alen","gaps","evalue"))
#
# [0.3] compute the cover and filter by thresholds and alignment
blast_DNA$qcover <- blast_DNA$alen/blast_DNA$qlen
keepers_DNA <- blast_DNA[blast_DNA$qcover < 0.2, ]$query
# > length(keepers_DNA)
# [1] 97860
#
no_hits_DNA <- setdiff(llist$V1, blast_DNA$query)
# > length(no_hits_DNA)
# [1] 209520
#
# > summary(blast_DNA[blast_DNA$query %in% keepers_DNA,]$qcover)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.002141 0.037550 0.069750 0.080170 0.115400 0.199900
#
# [0.4] do some set algebra to get the list of "interesting contigs"
keepers_set <- unique(c(keepers_DNA, no_hits_DNA))
# > length(keepers_set)
# [1] 307380
# 
# all these we need
write.table(keepers_set,"keepers_refseq_DNA.list",row.names = F,col.names = F,quote = F)
#
# include_mf 
#
# $ include_mf contigs_after_refseqPEP.fa contigs_after_refseqDNA.fa keepers_refseq_DNA.list
# phase 0 - getting buffers & file handlers ... 
# . input fasta: contigs_after_refseqPEP.fa 
# . output fasta: contigs_after_refseqDNA.fa 
# . include list: keepers_refseq_DNA.list 
# . skipped sequences list: keepers_refseq_DNA.list.skipped 
# phase 1 - counting sequences in the include list ... 
# ... 307380 sequences is about to be included, reading the list...
# ... 307380 sequences loaded, sorting ...
# phase 2 - filtering input ... 
# finished - included 307380 sequences
