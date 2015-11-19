#
# NOTE! THIS DELETES FOLDERS!
#
#
# set folders
drap_out_folder <-
  "/work2/project/sigenae/Project_RNA_Refinery/runDrapOut"
archive_dir <-
  "/work2/project/sigenae/Project_RNA_Refinery/assembled"
#
# set the drap out folder as the working dir
setwd(drap_out_folder)
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
dirs <- list.dirs(drap_out_folder, full.names = T)
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
      path, "*ASSEMBLY_COMPLETE*", full.names = FALSE, recursive = FALSE, ignore.case = TRUE
    )
  if (length(all[!file.info(all)$isdir]) == 1) {
    # is complete
    return(all)
  }else {
    # a candidate for a restart
    return("ASSEMBLY_INCOMPLETE")
  }
}
#
library(plyr)
#
assemblys$complete <-
  daply(assemblys, .(accession), function(x) {
    find_complete(x$fullpath)
  })
#
# filter out only complete assemblies
assemblys <-
  assemblys[grep("ASSEMBLY_COMPLETE", assemblys$complete),]
#
# populate links for fpkm1
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
#
assemblys$fpkm1 <-
  daply(assemblys, .(accession), function(x) {
    find_file(x$fullpath, "^transcripts_fpkm_1.fa$")
  })
assemblys$fpkm3 <-
  daply(assemblys, .(accession), function(x) {
    find_file(x$fullpath, "^transcripts_fpkm_3.fa$")
  })
assemblys$xprs <- daply(assemblys, .(accession), function(x) {
  find_file(x$fullpath, "^results.xprs$", recursive = TRUE)
})
assemblys <- assemblys[(length(assemblys$fpkm1) > 1) &
                         (length(assemblys$fpkm3) > 1) & (length(assemblys$xprs) > 1),]
str(assemblys)
#
#
global_res = rep("",7)
for (assembly_accession in assemblys$accession) {
  res <- c(assembly_accession)
  out_folder <- paste(archive_dir, assembly_accession, sep = "/")
  dir.create(out_folder, showWarnings = FALSE, recursive = TRUE)
  
  in_file <-
    assemblys[assemblys$accession == assembly_accession,]$fpkm1
  out_file <- paste(out_folder, basename(in_file), sep = "/")
  file.copy(in_file, out_folder)
  res <- cbind(res, out_file, file.info(out_file)$size)
  
  in_file <-
    assemblys[assemblys$accession == assembly_accession,]$fpkm3
  out_file <- paste(out_folder, basename(in_file), sep = "/")
  file.copy(in_file, out_folder)
  res <- cbind(res, out_file, file.info(out_file)$size)
  
  in_file <-
    assemblys[assemblys$accession == assembly_accession,]$xprs
  out_file <- paste(out_folder, basename(in_file), sep = "/")
  file.copy(in_file, out_folder)
  res <- cbind(res, out_file, file.info(out_file)$size)
  print(paste(res))
  global_res = rbind(global_res, res)
}
#
# ASSEMBLY FOLDERS DELETION
#
#rm_res <-
#  alply(assemblys[1,]$fullpath, 1, function(x) {
#    unlink(x, recursive = TRUE)
#  })