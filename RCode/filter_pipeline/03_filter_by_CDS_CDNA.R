require(stringr)
require(plyr)
require(dplyr)
require(data.table)
#
# [0.1] take care about the input file names
#
wfolder <- "/work2/project/sigenae/Project_RNA_Refinery/pipelineOut_embryo/"
raw_contigs_llist <- paste(wfolder, "raw_contigs.llist", sep = "")
blat_out_CDS <- paste(wfolder, 
      "blatOut_Galgal4_cds/raw_contigs_embryo.fa.Gallus_gallus.Galgal4.cds.all.fa.blat.best.tsv", sep = "")
blat_out_CDNA <- paste(wfolder, 
      "blatOut_Galgal4_cdna/raw_contigs_embryo.fa.Gallus_gallus.Galgal4.cdna.all.fa.blat.best.tsv", sep = "")
#
# [0.2] read the input
#
# makellist raw_contigs.fa >raw_contigs.llist
llist <- read.table(
    contigs_llist, header = FALSE, sep = " ", quote = "",
    stringsAsFactors = FALSE,comment.char = "",
    colClasses = c("character","integer"))
# > dim(llist)
# [1] 1442455       2
#
blat_CDS <- fread(blat_out_CDS)
# > dim(blat_CDS)
# [1] 860336     17
#
blat_CDNA <- fread(blat_out_CDNA)
# > dim(blat_CDNA)
# [1] 1009645      17
#
#
# [0.3] make lists of:
#       - those contigs whose best hit covers less than 20% of reference
#       - those contigs that have no hits
# keepers_CDS <- blat_CDS[blat_CDS$"%qCoverage" < 20.0,]
# keepers_CDNA <- blat_CDNA[blat_CDNA$"%qCoverage" < 20.0,]
#
keepers_CDS <- blat_CDS[as.numeric(gsub("%", "", blat_CDS$qCoverage)) < 20.0,]$qName
# > length(keepers_CDS)
# [1] 1448
#
keepers_CDNA <- blat_CDNA[as.numeric(gsub("%", "", blat_CDNA$qCoverage)) < 20.0,]$qName
# > length(keepers_CDNA)
# [1] 52140
#
no_hits_CDS <- setdiff(llist$V1, blat_CDS$qName)
no_hits_CDNA <- setdiff(llist$V1, blat_CDNA$qName)
# > length(no_hits_CDS)
# [1] 582119
# > length(no_hits_CDNA)
# [1] 432810
#
# [0.4] do some set algebra to get the list of "interesting contigs"
# 
# this is the clear winner -- no hits at all
no_hits_set <- intersect(no_hits_CDS, no_hits_CDNA)
# > length(no_hits_set)
# [1] 431286
# 
# this one is also a winning set -- hit is less 20% coverage in CDS and CDNA
keepers_set <- intersect(keepers_CDS, keepers_CDNA)
# > length(keepers_set)
# [1] 62759
#
# these two are good ones also
keepers_CDNA <- intersect(keepers_CDNA, no_hits_CDS)
keepers_CDS <- intersect(keepers_CDS, no_hits_CDNA)
# > length(keepers_CDNA)
# [1] 52140
# > length(keepers_CDS)
# [1] 1448
#
# final set
keepers <- unique( c(no_hits_set, keepers_set, keepers_CDNA, keepers_CDS))
# > length(keepers)
# [1] 547633
# 
# > summary( as.numeric(gsub("%", "", blat_CDS[blat_CDS$qName %in% keepers, ]$qCoverage ) ))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.300   4.400   8.200   8.863  12.900  19.900 
#
# > summary( as.numeric(gsub("%", "", blat_CDNA[blat_CDNA$qName %in% keepers, ]$qCoverage ) ))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.300   3.900   7.600   8.634  13.000  19.900
# 
# > length( which(blat_CDNA$qName  %in% keepers ))
# [1] 114899
# 
# > length( which(blat_CDS$qName  %in% keepers ))
# [1] 64207
# 
# > length( no_hits_set[no_hits_set %in% keepers] )
# [1] 431286
#
# [0.5] write down the list of "keepers"
# all these we need
write.table(keepers,"keepers_cds_cdna.list",row.names = F,col.names = F,quote = F)
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
