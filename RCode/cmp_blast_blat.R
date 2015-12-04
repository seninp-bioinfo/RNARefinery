require(stringr)
require(reshape)
require(plyr)
require(dplyr)
require(data.table)
#
require(ggplot2)
require(gridExtra)
require(scales)
require(Cairo)
#
blat = fread(input = "~/tmp/refinery/blat_cdna_refseq_pep.best.tsv")
setnames(blat, gsub("%","",names(blat)))
blast = fread(input = "~/tmp/refinery/contigs_after_refseq_pep.best.tsv")
#
str(blat)
dd = full_join(blat, blast, by=c("tName"="V1"))
length(which(is.na(dd$qName)))/9658
length(which(is.na(dd$V1)))

dd = left_join(as.data.frame(blat), as.data.frame(blast), by = c("tName" = "V1"))
length(which(dd$qName != dd$V2))
length(dd$tName)
2128/length(dd$tName)
ggplot(data=blast, aes(x=V4)) + geom_density()
ggplot(data=blat, aes(x=identity)) + geom_density()

par(mfrow=c(1,2))
hist(blast$V4)
hist(blat$identity)

head(blat$tName)
head(blast$V1)
