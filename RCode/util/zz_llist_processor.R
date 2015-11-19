require(RMySQL)
#
require(rentrez)
require(XML)
#
require(reshape)
require(plyr)
require(dplyr)
#
require(stringr)
#
require(ggplot2)
require(grid)
require(gridExtra)
require(scales)
require(Cairo)
#
properties = read.table("/home/psenin/.rnarefinery/db.properties", header = F, as.is = T)
db_credentials = data.frame(t(properties$V2), stringsAsFactors = F); names(db_credentials) = properties$V1
session <- dbConnect(MySQL(), host = db_credentials$host, dbname = db_credentials$db, 
                     user = db_credentials$user, password = db_credentials$password)
#
# raw contigs
#
contigs_summary <- dbGetQuery(session, 
  paste("select sr.run_accession, sr.tissue, count(rc.name) as raw_contigs from raw_contigs rc
  join sra_runs sr on rc.run_id=sr.id group by rc.run_id"))
#
llist_summary <- function(llist_name) {
  llist <- read.table(llist_name, header = FALSE, sep = " ", quote = "",
    stringsAsFactors = FALSE,comment.char = "", colClasses = c("character","integer"))
  llist$run_accession <- str_extract(llist$V1, "[^_]*")
  sum <- ddply(llist,.(run_accession), summarize, contigs = length(run_accession))
  sum
}  
#
# RAW
#
llist_dd1 = llist_summary("/work2/project/sigenae/Project_RNA_Refinery/pipelineOut/raw_contigs.llist")
llist_dd2 = llist_summary("/work2/project/sigenae/Project_RNA_Refinery/pipelineOut_embryo/raw_contigs.llist")
dd = rbind(llist_dd1, llist_dd2)
dm = merge(contigs_summary, dd)
#
# CDS and CDNA
#
llist_dd1 = llist_summary("/work2/project/sigenae/Project_RNA_Refinery/pipelineOut/contigs_after_CDS_cDNA.llist")
llist_dd2 = llist_summary("/work2/project/sigenae/Project_RNA_Refinery/pipelineOut_embryo/contigs_after_CDS_cDNA.llist")
dd = rbind(llist_dd1, llist_dd2)
names(dd) <- c("run_accession", "contigs_after_cds_cdna")
dm2 = merge(dm, dd)
#
# REFSEQ PEP
#
llist_dd1 = llist_summary("/work2/project/sigenae/Project_RNA_Refinery/pipelineOut/contigs_after_PEP.llist")
llist_dd2 = llist_summary("/work2/project/sigenae/Project_RNA_Refinery/pipelineOut_embryo/contigs_after_refseqPEP.llist")
dd = rbind(llist_dd1, llist_dd2)
names(dd) <- c("run_accession", "contigs_after_refseq_pep")
dm3 = merge(dm2, dd)
#
# REFSEQ PEP
#
llist_dd1 = llist_summary("/work2/project/sigenae/Project_RNA_Refinery/pipelineOut/contigs_after_refseqDNA.llist")
llist_dd2 = llist_summary("/work2/project/sigenae/Project_RNA_Refinery/pipelineOut_embryo/contigs_after_refseqDNA.llist")
dd = rbind(llist_dd1, llist_dd2)
names(dd) <- c("run_accession", "contigs_after_refseq_dna")
dm4 = merge(dm3, dd)
#
#
dm4$tissue = as.factor(dm4$tissue)
dm4$run_accession = as.factor(dm4$run_accession)
#
dm_melt <- melt(dm4, id.var=c("run_accession", "tissue"))
p_tissue <- ggplot(data=dm_melt, aes(x=variable, y=value, color=tissue, group=run_accession)) + geom_line() +
  theme_bw() + scale_y_continuous("Contigs") + ggtitle("Filtering process by tissue")

p_accession <- ggplot(data=dm_melt, aes(x=variable, y=value, color=run_accession, group=run_accession)) + 
  geom_line() +  theme_bw() + scale_y_continuous("Contigs") + ggtitle("Filtering process by run accession") +
  guides(color=guide_legend(ncol=3))

CairoPDF(width = 16, height = 24, file = "filtering.pdf", bg = "transparent")
print(grid.arrange(p_tissue, p_accession, ncol=1))
dev.off()

