require(data.table)
#
options(echo = TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(paste("arguments supplied:",args))
#
# [0.1] take care about the input file names
contigs_llist <- args[1]
best_tsv <- args[2]
keepers_list <- args[3]
#
# [0.2] read the input
#
llist <- read.table(contigs_llist, header = FALSE, sep = " ", quote = "",  stringsAsFactors = FALSE,
                    comment.char = "", colClasses = c("character","integer"))
tsv <- fread(best_tsv)
#
# [0.3] filter the input by coverage
#
cov_keepers <- tsv[tsv$"%qCoverage" < 20.0,]$qName
#
# [0.4] filter the input by no hits
#
no_hits_set <- setdiff(contigs_llist$V1, tsv$qName)
#
# [0.5] compose the set and save
#
keepers <- unique( c(no_hits_set, cov_keepers) )
write.table(keepers,keepers_list,row.names = F,col.names = F,quote = F)