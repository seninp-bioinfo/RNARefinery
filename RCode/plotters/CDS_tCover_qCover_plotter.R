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
# make a dataframe of columns:
#
dirs <- list.dirs(archive_dir, full.names = T)
assemblys <-  data.frame(accession = basename(dirs), fullpath = dirs, stringsAsFactors = FALSE)
#
# figure which of these are complete by checking if the desired filename is found
# this adds a column complete / incomplete to the data frame
#
find_complete <- function(path = ".", file_name = "blat_raw_cds.best.tsv") {
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
# filter out only complete assemblies
#
assemblys <- assemblys[grep("^complete$", assemblys$complete),]
#
# populate other data, for this define a find file function
#
find_file <- function(path = ".", filename = "", recursive = FALSE) {
  print(paste("searching ", filename, " in ", path))
  all <- list.files(
    path, filename, full.names = TRUE, recursive = recursive, ignore.case = TRUE )
  all <- all[!file.info(all)$isdir]
  if (length(all) > 0) {
    return(all)
  }else {
    return(NA)
  }
}
#
# get the raw contigs
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
##################################################################################
##################################################################################
#
contigs_nums = ldply(assemblys_contigs_llist, function(x){length(x$contig_name)})
setnames(contigs_nums,c("accession","contigs_num"))
assemblys = merge(assemblys,contigs_nums)
#
# what we want is a list of files which are in passed to the next step
#
assemblys$filter_cds_contigs_llist <- 
  daply(assemblys, .(accession), function(x){find_file(x$fullpath, "contigs_after_cds.llist", recursive=T)})
assemblys_filter_cds_llist <- dlply(assemblys, .(accession), function(x){
  print(paste(x$filter_cds_contigs_llist))
  list=fread(input=as.character(x$filter_cds_contigs_llist), header=F)
  setnames(list, c("contig_name", "length"))
  list
})
#
#
#
idx <- sample(1:301,30)
dd = assemblys[idx,]
dd = arrange(dd,tissue)
#
#
#
plots_tqcover = list()
for(i in 1:length(dd$accession)){
  acc = dd$accession[i]
  incontigs_llist = assemblys_contigs_llist[[acc]]
  outcontigs_llist = assemblys_filter_cds_llist[[acc]]
  tsv = assemblys_filter_cds[[acc]]
  
  exclude_list <- tsv[(tsv$"%identity" >= 75.0) & 
                        ((tsv$"%tCoverage" >= 75.0) | (tsv$"%qCoverage" >= 75.0)),]$qName
  cov_keepers <- setdiff(tsv$qName, exclude_list)
  
  no_hits_set <- setdiff(incontigs_llist$contig_name, tsv$qName)
  
  print(paste("keepers homology:",length(cov_keepers),"keepers no hits:",length(no_hits_set)))
  keepers <- unique( c(no_hits_set, cov_keepers) )
  
  setdiff(outcontigs_llist$contig_name, keepers)
}

Cairo(width = 1600, height = 1600, 
      file="after_75cutoff_tcoverage.pdf", type="pdf", pointsize=24, 
      bg = "transparent", canvas = "white", units = "px", dpi = 55)
do.call(grid.arrange,  plots_tqcover)
dev.off()


