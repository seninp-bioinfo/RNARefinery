require(reshape)
require(plyr)
require(dplyr)
require(stringr)
#
require(RSQLite)
require(SRAdb)
require(RMySQL)
#
require(ggplot2)
require(gridExtra)
require(scales)
require(Cairo)
#
# set folders
#
archive_dir <- "/work2/project/sigenae/Project_RNA_Refinery/assembled"
#
# set the working dir
setwd(archive_dir)
#
# list all the folders within the working directory
list.dirs <-
  function(path = ".", pattern = NULL, all.dirs = FALSE, full.names = FALSE, ignore.case = FALSE) {
    # use full.names=TRUE to pass to file.info
    all <-
      list.files(path, pattern, all.dirs, full.names = TRUE, recursive = FALSE, ignore.case)
    dirs <- all[file.info(all)$isdir]
    # determine whether to return full names or just dir names
    if (isTRUE(full.names))
      return(dirs)
    else
      return(basename(dirs))
  }
#
# make a dataframe of an SRA accession and the full path
dirs <- list.dirs(archive_dir, full.names = T)
assemblys <-
  data.frame(
    accession = basename(dirs), fullpath = dirs, stringsAsFactors = FALSE
  )
#
# figure which of these are complete
find_complete <- function(path = ".") {
  print(paste("searching ", path))
  all <-
    list.files(
      path, "contigs_after_ncrna.fa", full.names = TRUE, recursive = FALSE, ignore.case = TRUE
    )
  if (length(all[!file.info(all)$isdir]) == 1) {
    # is complete
    return(all)
  }else {
    # a candidate for a restart
    return("INCOMPLETE")
  }
}
#
assemblys$complete <-
  daply(assemblys, .(accession), function(x) {
    find_complete(paste(x$fullpath,"/filters",sep = ""))
  })
#
# filter out only complete assemblies
assemblys <-
  assemblys[grep("contigs_after_ncrna", assemblys$complete),]
#
# get connected to local DB copy
sqlfile <- "/work2/project/sigenae/Project_RNA_Refinery/SRAmetadb.sqlite"
sra_con <- dbConnect(SQLite(),sqlfile)
# select RUNS
sra_info <- function(run_accession){
  rs <- dbGetQuery(sra_con,paste(
    "SELECT sra.* FROM sra WHERE run_accession=",shQuote(run_accession),sep = ""))
  rs
}  
#
# retrieve SPOTS from SRA
assemblys$spots <- daply(assemblys, .(accession), function(x){sra_info(x$accession)$spots})
#
# populate other data
find_file <- function(path = ".", filename = "", recursive = FALSE) {
  print(paste("searching ", filename, " in ", path))
  all <-
    list.files(
      path, filename, full.names = TRUE, recursive = recursive, ignore.case = TRUE
    )
  all <- all[!file.info(all)$isdir]
  if (length(all) > 0) {
    # is complete
    return(all)
  }else {
    # a candidate for a restart
    return("")
  }
}
#
# retrieve raw contigs data
assemblys$contigs_llist <- daply(assemblys, .(accession), function(x){
  find_file(x$fullpath, "raw_contigs.llist", recursive=T)})
filter_raw_llist <- dlply(assemblys, .(accession), function(x){
  llist <- read.table(x$contigs_llist, header = FALSE, sep = " ", quote = "",  stringsAsFactors = FALSE,
                      comment.char = "", colClasses = c("character","integer"))
  llist
})
assemblys_contigs = ldply(filter_raw_llist,function(x){length(x$V1)})
names(assemblys_contigs) <- c("accession","raw_contigs")
assemblys = merge(assemblys,assemblys_contigs)
#
# get the TISSUE from mysql
#
#
properties = read.table("/home/psenin/.rnarefinery/db.properties", header = F, as.is = T)
db_credentials = data.frame(t(properties$V2), stringsAsFactors = F); names(db_credentials) = properties$V1
session <- dbConnect(MySQL(), host = db_credentials$host, dbname = db_credentials$db, 
                     user = db_credentials$user, password = db_credentials$password)
a_tissue <- ddply(assemblys, .(accession), function(x){
  tissue = dbGetQuery(session, paste("select tissue from sra_runs where run_accession=",
      shQuote(x$accession)))
  tissue
})
names(a_tissue) <- c("accession","tissue")
assemblys = merge(assemblys,a_tissue)
#
#
tissue_table=as.data.frame(table(assemblys$tissue))
names(tissue_table)=c("tissue","frequency")
tissue_table$str = paste(tissue_table$tissue," (",tissue_table$frequency,")",sep="")
tissue_label = ddply(assemblys,.(accession),function(x){tissue_table[tissue_table$tissue==x$tissue,]$str})
names(tissue_label) = c("accession","tissue_label")
assemblys = merge(assemblys, tissue_label)
p=ggplot(assemblys,aes(x=as.factor(tissue_label),y=raw_contigs))+geom_boxplot(aes(fill=as.factor(tissue_label))) +
  guides(fill=guide_legend(ncol=2)) + theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
  ggtitle("Assembled contigs variance per tissue")
Cairo(width = 1600, height = 900, 
      file="tissue_raw_contigs.pdf", 
      type="pdf", pointsize=12, 
      bg = "transparent", canvas = "white", units = "px", dpi = 82)
print(p)
dev.off()