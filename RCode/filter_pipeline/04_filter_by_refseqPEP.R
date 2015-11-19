require(stringr)
require(plyr)
require(dplyr)
require(data.table)
#
# [0.1] take care about the input file names
#
wfolder <- "/work2/project/sigenae/Project_RNA_Refinery/pipelineOut_embryo/"
contigs_llist <- paste(wfolder, "contigs_after_CDS_cDNA.llist", sep = "")
blastx_PEP <- paste(wfolder, "contigs_after_CDS_CDNA.refseqPEP.blasted.table", sep = "")
#
# [0.2] read the input
#
llist <- read.table(
    contigs_llist, header = FALSE, sep = " ", quote = "",
    stringsAsFactors = FALSE,comment.char = "",
    colClasses = c("character","integer"))
# > dim(llist)
# [1] 547633      2
# 
blast_PEP <- fread(blastx_PEP)
# > dim(blast_PEP)
# [1] 193201      9
#
setnames(blast_PEP, c("query","hit","desc","identity","bits","qlen","alen","gaps","evalue"))
#
# [0.3] compute the cover and filter by thresholds and alignment
blast_PEP$qcover <- blast_PEP$alen/blast_PEP$qlen
keepers_PEP <- blast_PEP[blast_PEP$qcover < 0.2, ]$query
# > length(keepers_PEP)
# [1] 119488
#
no_hits_PEP <- setdiff(llist$V1, blast_PEP$query)
# > length(no_hits_PEP)
# [1] 354432
#
# > summary(blast_PEP[blast_PEP$query %in% keepers_PEP,]$qcover)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.002011 0.043080 0.083170 0.092430 0.141300 0.200000
#
# [0.4] do some set algebra to get the list of "interesting contigs"
keepers_set <- unique(c(keepers_PEP, no_hits_PEP))
# > length(keepers_set)
# [1] 473920
# 
# all these we need
write.table(keepers_set,"keepers_refseq_PEP.list",row.names = F,col.names = F,quote = F)
#
# include_mf 
#
# $ include_mf contigs_after_CDS_cDNA.fa contigs_after_refseqPEP.fa keepers_refseq_PEP.list
# phase 0 - getting buffers & file handlers ... 
# . input fasta: contigs_after_CDS_cDNA.fa 
# . output fasta: contigs_after_refseqPEP.fa 
# . include list: keepers_refseq_PEP.list 
# . skipped sequences list: keepers_refseq_PEP.list.skipped 
# phase 1 - counting sequences in the include list ... 
# ... 473920 sequences is about to be included, reading the list...
# ... 473920 sequences loaded, sorting ...
# phase 2 - filtering input ... 
# finished - included 473920 sequences