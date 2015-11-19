require(stringr)
require(plyr)
require(dplyr)
require(data.table)
#
# [0.1] take care about the input file names
#
wfolder <- "/work2/project/sigenae/Project_RNA_Refinery/pipelineOut_embryo/"
raw_contigs_llist <- paste(wfolder, "contigs_after_refseqDNA.llist", sep = "")
blat_out_genomic <- paste(wfolder, 
      "blatOut_Galgal4_genomic/contigs_after_refseqDNA.fa.Gallus_gallus.Galgal4.dna_sm.toplevel.fa.32915.blat.best.tsv", sep = "")
blat_out_ncRNA <- paste(wfolder, 
      "blatOut_Galgal4_ncRNA/contigs_after_refseqDNA.fa.Gallus_gallus.Galgal4.ncrna.fa.33087.blat.best.tsv", sep = "")
#
# [0.2] read the input
#
# makellist raw_contigs.fa >raw_contigs.llist
llist <- read.table(raw_contigs_llist, header = FALSE, sep = " ", quote = "",
    stringsAsFactors = FALSE,comment.char = "", colClasses = c("character","integer"))
# > dim(llist)
# [1] 307380      2
#
blat_genomic <- fread(blat_out_genomic)
# > dim(blat_genomic)
# [1] 266773     19
#
blat_ncRNA <- fread(blat_out_ncRNA)
# > dim(blat_ncRNA)
# [1] 4940   19
#
#
# [0.3] make lists of:
#       - those contigs whose best hit covers less than 20% of reference
#       - those contigs that have no hits
keepers_genomic <- blat_genomic[blat_genomic$"%qCoverage" < 20.0,]$qName
# > length(keepers_genomic)
# [1] 19255
keepers_ncRNA <- blat_ncRNA[blat_ncRNA$"%qCoverage" < 20.0,]$qName
# > length(keepers_ncRNA)
# [1] 4164
#
no_hits_genomic <- setdiff(llist$V1, blat_genomic$qName)
no_hits_ncRNA <- setdiff(llist$V1, blat_ncRNA$qName)
# > length(no_hits_genomic)
# [1] 40607
# > length(no_hits_ncRNA)
# [1] 302440
#
# [0.4] do some set algebra to get the list of "interesting contigs"
# 
# this is the clear winner -- no hits at all
no_hits_set <- intersect(no_hits_genomic, no_hits_ncRNA)
# > length(no_hits_set)
# [1] 40547
# 
# this one is also a winning set -- hit is less 20% coverage in CDS and CDNA
keepers_set <- intersect(keepers_genomic, keepers_ncRNA)
# > length(keepers_set)
# [1] 134
#
# these two are good ones also
keepers_set_genomic <- intersect(keepers_genomic, no_hits_ncRNA)
keepers_set_ncRNA <- intersect(keepers_ncRNA, no_hits_genomic)
# > length(keepers_genomic)
# [1] 19255
# > length(keepers_ncRNA)
# [1] 4164
#
# final set
keepers <- unique( c(no_hits_set, keepers_set, keepers_set_genomic, keepers_set_ncRNA))
# > length(keepers)
# [1] 59862
# 
# > summary( blat_genomic[blat_genomic$qName %in% keepers, ]$"%qCoverage" )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.86    6.26    9.84   10.29   14.22   19.99 

# > summary( blat_ncRNA[blat_ncRNA$qName %in% keepers, ]$"%qCoverage" )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.040   8.208   9.640  10.020  11.960  19.030# 
#
# [0.5] write down the list of "keepers"
# all these we need
write.table(setdiff(llist$V1, keepers),"keepers_genomic.list",row.names = F,col.names = F,quote = F)
#
# include_mf 
#
# $ include_mf raw_contigs_embryo.fa contigs_after_CDS_cDNA.fa keepers_cds_cdna.list
# phase 0 - getting buffers & file handlers ... 
# . input fasta: raw_contigs_embryo.fa 
# . output fasta: contigs_after_CDS_cDNA.fa 
# . include list: keepers_cds_cdna.list 
# . skipped sequences list: keepers_cds_cdna.list.skipped 
# phase 1 - counting sequences in the include list ... 
# ... 547633 sequences is about to be included, reading the list...
# ... 547633 sequences loaded, sorting ...
# phase 2 - filtering input ... 
# finished - included 547633 sequences
