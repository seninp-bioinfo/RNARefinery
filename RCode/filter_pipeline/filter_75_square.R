require(data.table)
#
options(echo = TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(paste("arguments supplied:",args))
#
# [0.1] take care about the input file names
llist_fname <- args[1]
tsv_fname <- args[2]
keepers_fname <- args[3]
#
# [0.2] read the input
#
llist <- read.table(llist_fname, header = FALSE, sep = " ", quote = "",  stringsAsFactors = FALSE,
                    comment.char = "", colClasses = c("character","integer"))
tsv <- fread(tsv_fname)
#
# [0.3] filter the input by coverage
#
exclude_list <- tsv[(tsv$"%identity" >= 75.0) & 
      ((tsv$"%tCoverage" >= 75.0) | (tsv$"%qCoverage" >= 75.0)),]$qName
cov_keepers <- setdiff(tsv$qName, exclude_list)
#summary(tsv[tsv$qName %in% cov_keepers, ]$"%qCoverage" )
#summary(tsv[tsv$qName %in% cov_keepers, ]$"%tCoverage" )
#summary(tsv[tsv$qName %in% cov_keepers, ]$"%identity" )
#
# [0.4] filter the input by no hits
#
no_hits_set <- setdiff(llist$V1, tsv$qName)
#
# [0.5] compose the set and save
#
print(paste("keepers homology:",length(cov_keepers),"keepers no hits:",length(no_hits_set)))
keepers <- unique( c(no_hits_set, cov_keepers) )
#
write.table(keepers,keepers_fname,row.names = F,col.names = F,quote = F)
