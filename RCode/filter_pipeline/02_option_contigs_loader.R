# set folders
archive_dir <- "/work2/project/sigenae/Project_RNA_Refinery/assembled"
#
# set the drap out folder as the working dir
setwd(archive_dir)
#
# list all the folders within the working directory
list.dirs <- function(path = ".", pattern = NULL, all.dirs = FALSE, full.names = FALSE, ignore.case = FALSE) {
  # use full.names=TRUE to pass to file.info
  all <- list.files(path, pattern, all.dirs, full.names = TRUE, recursive = FALSE, ignore.case)
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
#
#
library(seqinr)
require(RMySQL)
properties = read.table("/home/psenin/.rnarefinery/db.properties", header = F, as.is = T)
db_credentials = data.frame(t(properties$V2), stringsAsFactors = F); names(db_credentials) = properties$V1
session <- dbConnect(MySQL(), host = db_credentials$host, dbname = db_credentials$db, 
                     user = db_credentials$user, password = db_credentials$password)
#
for (dir in dirs) {
  accession <- basename(dir)
  print(paste("** processing", accession))
  
  sra_run <- dbGetQuery(session, paste("select * from sra_runs where run_accession=",shQuote(accession)))
  fin_name <- paste(dir,"transcripts_fpkm_3.fa",sep = "/")
  fin <- file.info(fin_name)
  
  if (is.na(fin$size) || fin$size < 1 || is.na(sra_run$id)) {
    stop(paste("Error with the run ", accession, dir))
  }
  sequences <- read.fasta(fin_name, as.string = T)  
  #
  print(paste("  .. read", length(names(sequences)), "sequences"))
  #
  ctr <- 0
  for (seqname in names(sequences)) {
    
    dbGetQuery(session, paste("INSERT IGNORE INTO raw_contigs (run_id, name, sequence) VALUES (",
                              shQuote(sra_run$id), ",",
                              shQuote(paste(accession,seqname, sep = "_")),",",
                              shQuote(sequences[[seqname]][1]),")"))
    ctr <- ctr + 1
    if (ctr %% 1000 == 0) {
      print(paste("  .. loaded ", ctr, " sequences for the run ", accession, ", id ", sra_run$id))
    }
  }
}

