# table_maker
#
require(stringr)
require(reshape)
require(plyr)
require(dplyr)
require(data.table)
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
setwd(archive_dir)
#
# list all the folders within the working directory
#
list.dirs <-
  function(path = ".", pattern = NULL, all.dirs = FALSE, full.names = FALSE, ignore.case = FALSE) {
    # use full.names=TRUE to pass to file.info
    all <- list.files(path, pattern, all.dirs, full.names = TRUE, recursive = FALSE, ignore.case)
    dirs <- all[file.info(all)$isdir]
    if (isTRUE(full.names))
      return(dirs)
    else
      return(basename(dirs))
  }
#
""
# make a dataframe of columns:
#    1) SRA accession
#
dirs <- list.dirs(archive_dir, full.names = T)
assemblys <-  data.frame(accession = basename(dirs), fullpath = dirs, stringsAsFactors = FALSE)
#
# figure which of these are complete by checking if the desired filename is found
# this adds a column complete / incomplete to the data frame
#
find_complete <- function(path = ".", file_name = "contigs_after_ncrna.fa") {
  print(paste("searching ", path, "for", file_name))
  all <-
    list.files(
      path, file_name, full.names = TRUE, recursive = FALSE, ignore.case = TRUE
    )
  if (length(all[!file.info(all)$isdir]) == 1) {
    # is complete
    return("complete")
  }else {
    # a candidate for a restart
    return("incomplete")
  }
}
#
assemblys$complete <- daply(assemblys, .(accession), function(x) {
  find_complete(paste(x$fullpath,"/filters",sep = ""))
})
#
#
# attach the run tissue description
#
properties = read.table("/home/psenin/.rnarefinery/db.properties", header = F, as.is = T)
db_credentials = data.frame(t(properties$V2), stringsAsFactors = F); names(db_credentials) = properties$V1
session <- dbConnect(MySQL(), host = db_credentials$host, dbname = db_credentials$db, 
                     user = db_credentials$user, password = db_credentials$password)
a_tissue <- ddply(assemblys, .(accession), function(x){
  tissue = dbGetQuery(session, paste("select tissue from sra_runs where run_accession=", shQuote(x$accession)))
  tissue
})
names(a_tissue) <- c("accession","tissue")
assemblys = merge(assemblys,a_tissue)
#
#
# filter out only complete assemblies
assemblys <- assemblys[grep("^complete$", assemblys$complete),]
#
# populate other data, for this define a find file function
#
find_file <- function(path = ".", filename = "", recursive = FALSE) {
  print(paste("searching ", filename, " in ", path))
  all <-list.files(
    path, filename, full.names = TRUE, recursive = recursive, ignore.case = TRUE )
  all <- all[!file.info(all)$isdir]
  if (length(all) > 0) {
    return(all)
  }else {
    return(NA)
  }
}
#
# get the contigs 
#
assemblys$contigs_llist <- daply(assemblys, .(accession), function(x){find_file(x$fullpath, "raw_contigs.llist", recursive=T)})
#
# get the contigs and the first CDS alignment
#
assemblys$filter0_cds_tsv <- daply(assemblys, .(accession), function(x){find_file(x$fullpath, "blat_raw_cds.best.tsv", recursive=T)})
#
# get the CDS tsv and the LLIST
#
assemblys_filter_cds <- dlply(assemblys, .(accession), function(x){
  print(paste(x$filter0_cds_tsv))
  tsv = fread(input = as.character(x$filter0_cds_tsv))
  tsv
})
#
assemblys_contigs_llist <- dlply(assemblys, .(accession), function(x){
  print(paste(x$contigs_llist))
  list=fread(input=as.character(x$contigs_llist), header=F)
  setnames(list, c("contig_name", "length"))
  list
})
#
#
contigs_nums = ldply(assemblys_contigs_llist, function(x){length(x$contig_name)})
setnames(contigs_nums,c("accession","contigs_num"))
assemblys = merge(assemblys,contigs_nums)
#
#
contigs_over75 = ldply(assemblys_filter_cds, function(x){
  setnames(x, gsub("%","",names(x)))
  keepers = x[x$tCoverage<75 & x$identity<75,]
  length(keepers$qName)
})
setnames(contigs_over75,c("accession","keepers_cds_75"))
assemblys = merge(assemblys,contigs_over75)

write.table(assemblys,"75cutoff_experiment.csv",row.names=F,col.names=T,quote="\"",se)