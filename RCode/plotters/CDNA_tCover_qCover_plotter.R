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
archive_dir <-
  "/work2/project/sigenae/Project_RNA_Refinery/assembled"
setwd(archive_dir)
#
# list all the folders within the working directory
#
list.dirs <-
  function(path = ".", pattern = NULL, all.dirs = FALSE, full.names = FALSE, ignore.case = FALSE) {
    # use full.names=TRUE to pass to file.info
    all <-
      list.files(path, pattern, all.dirs, full.names = TRUE, recursive = FALSE, ignore.case)
    dirs <- all[file.info(all)$isdir]
    if (isTRUE(full.names))
      return(dirs)
    else
      return(basename(dirs))
  }
#
# make a dataframe of columns:
#    1) SRA accession
#
dirs <- list.dirs(archive_dir, full.names = T)
assemblys <-
  data.frame(
    accession = basename(dirs), fullpath = dirs, stringsAsFactors = FALSE
  )
#
# figure which of these are complete by checking if the desired filename is found
# this adds a column complete / incomplete to the data frame
#
find_complete <-
  function(path = ".", file_name = "contigs_after_ncrna.fa") {
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
  find_complete(paste(x$fullpath,"/filters/",sep = "") ,"contigs_after_cdna.llist")
})
#
# attach the tissue description
#
properties = read.table("/home/psenin/.rnarefinery/db.properties", header = F, as.is = T)
db_credentials = data.frame(t(properties$V2), stringsAsFactors = F); names(db_credentials) = properties$V1
session <-
  dbConnect(MySQL(), host = db_credentials$host, dbname = db_credentials$db,
    user = db_credentials$user, password = db_credentials$password)
  a_tissue <- ddply(assemblys, .(accession), function(x) {
  tissue = dbGetQuery(
    session, paste("select tissue from sra_runs where run_accession=", shQuote(x$accession)))
  tissue
})
names(a_tissue) <- c("accession","tissue")
assemblys = merge(assemblys,a_tissue)
#
# filter out only complete assemblies
#
assemblys <- assemblys[grep("^complete$", assemblys$complete),]
#
# get connected to local DB copy and get some SRA info
#
sqlfile <-
  "/work2/project/sigenae/Project_RNA_Refinery/SRAmetadb.sqlite"
sra_con <- dbConnect(SQLite(),sqlfile)
#
# define the SRA query function
#
sra_info <- function(run_accession) {
  rs <-
    dbGetQuery(sra_con,paste(
      "SELECT sra.* FROM sra WHERE run_accession=",shQuote(run_accession),sep = ""
    ))
  rs
}
assemblys$reads <-
  daply(assemblys, .(accession), function(x) {
    sra_info(x$accession)$spots
  })
#
# populate other data, for this define a find file function
#
find_file <-
  function(path = ".", filename = "", recursive = FALSE) {
    print(paste("searching ", filename, " in ", path))
    all <- list.files(
      path, filename, full.names = TRUE, recursive = recursive, ignore.case = TRUE
    )
    all <- all[!file.info(all)$isdir]
    if (length(all) > 0) {
      return(all)
    }else {
      return(NA)
    }
  }
#
# Get the raw contigs LLIST
#
assemblys$contigs_llist <-
  daply(assemblys, .(accession), function(x) {
    find_file(x$fullpath, "raw_contigs.llist", recursive = T)
  })
assemblys_contigs_llist <- dlply(assemblys, .(accession), function(x){
  print(paste(x$contigs_llist))
  list = fread(input = as.character(x$contigs_llist), header = F)
  setnames(list, c("contig_name", "length"))
  list
})
#
##################################################################################
####### A SAMPLE FOR PLOTTING DEFINED ############################################
##################################################################################
#
idx <- sample(1:301,16)
idx <- c(51,34,116,72,97,292,2,193,170,210,36,255,20,270,164,267)
dd = assemblys[idx,]
dd = arrange(dd,tissue)
#
##################################################################################
####### CDS ######################################################################
##################################################################################
#
# get the first CDS alignment
#
assemblys$filter0_cds_tsv <-
  daply(assemblys, .(accession), function(x) {
    find_file(x$fullpath, "blat_raw_cds.best.tsv", recursive = T)
  })
assemblys_filter_cds <- dlply(assemblys, .(accession), function(x){
  print(paste(x$filter0_cds_tsv))
  tsv = fread(input = as.character(x$filter0_cds_tsv))
  setnames(tsv, gsub("%","",names(tsv)))
  tsv
})
#
# plot the CDS alignment
#
plots_tqcover_cds_before = list()
for (i in 1:length(dd$accession)) {
  acc = dd$accession[i]
  #
  incontigs_llist = assemblys_contigs_llist[[acc]]
  tsv = assemblys_filter_cds[[acc]]
  no_hits = setdiff(incontigs_llist$contig_name, tsv$qName)
  title = paste("CDS:",dd$accession[i],dd$tissue[i],"\n",length(no_hits),"out of",
    length(incontigs_llist$contig_name), "contigs were not aligned")
  #
  plots_tqcover_cds_before[[i]] = ggplot(data = tsv, aes(x = tCoverage, y = qCoverage, 
    colour = identity)) +geom_jitter(alpha = 0.5) + geom_density2d() + theme(legend.position = "bottom") +
    ggtitle(title) + scale_x_continuous(limits=c(0,100)) + scale_y_continuous(limits=c(0,100)) +
    scale_colour_gradientn(
      name = "Identity:  ",limits = c(50,100),
      colours = c("red","yellow","green","lightblue","darkblue"),
      breaks = c(50,75,100),labels = c("low(50)","medium(75)","high(100)"),
      guide = guide_colorbar(
        title.theme = element_text(size = 14, angle = 0),title.vjust = 1,
        barheight = 0.6, barwidth = 6, label.theme = element_text(size = 10, angle = 0)
      )
    )
}
Cairo(
  width = 1200, height = 1200,
  file = "cds_alignment.pdf", type = "pdf", pointsize = 20,
  bg = "transparent", canvas = "white", units = "px", dpi = 60
)
do.call(grid.arrange,  plots_tqcover_cds_before)
dev.off()
#
# Get the out filter contigs LLIST
#
assemblys$contigs_after_cds_llist <-
  daply(assemblys, .(accession), function(x) {
    find_file(x$fullpath, "contigs_after_cds.llist", recursive = T)
  })
assemblys_contigs_after_cds_llist <- dlply(assemblys, .(accession), function(x){
  print(paste(x$contigs_after_cds_llist))
  list = fread(input = as.character(x$contigs_after_cds_llist), header = F)
  setnames(list, c("contig_name", "length"))
  list
})
#
# plot the out CDS filter figure
#
plots_tqcover_cds_after = list()
for (i in 1:length(dd$accession)) {
  acc = dd$accession[i]
  #
  incontigs_llist = assemblys_contigs_llist[[acc]]
  tsv = assemblys_filter_cds[[acc]]
  no_hits = setdiff(incontigs_llist$contig_name, tsv$qName)
  outcontigs_llist = assemblys_contigs_after_cds_llist[[acc]]
  tsv = filter(assemblys_filter_cds[[acc]], qName %in% outcontigs_llist$contig_name)
  title = paste("filtered CDS:",dd$accession[i],dd$tissue[i],"\n",length(no_hits),"out of",
    length(incontigs_llist$contig_name), "contigs were not aligned")
  #
  plots_tqcover_cds_after[[i]] = ggplot(data = tsv, aes(x = tCoverage, y = qCoverage, 
      colour = identity)) + geom_jitter(alpha = 0.5) + geom_density2d() + theme(legend.position = "bottom") +
    ggtitle(title) + scale_x_continuous(limits=c(0,100)) + scale_y_continuous(limits=c(0,100)) +
    scale_colour_gradientn(
      name = "Identity:  ",limits = c(50,100),
      colours = c("red","yellow","green","lightblue","darkblue"),
      breaks = c(50,75,100),labels = c("low(50)","medium(75)","high(100)"),
      guide = guide_colorbar(
        title.theme = element_text(size = 14, angle = 0),title.vjust = 1,
        barheight = 0.6, barwidth = 6, label.theme = element_text(size = 10, angle = 0)
      )
    )
}
Cairo(width = 1200, height = 1200, file = "cds_alignment_filtered.pdf", type = "pdf", pointsize = 20,
  bg = "transparent", canvas = "white", units = "px", dpi = 60 )
do.call(grid.arrange,  plots_tqcover_cds_after)
dev.off()
#
##################################################################################
####### cDNA #####################################################################
##################################################################################
#
assemblys$filter1_cdna_tsv <-
  daply(assemblys, .(accession), function(x) {
    find_file(x$fullpath, "blat_cds_cdna.best.tsv", recursive = T)
  })
assemblys_filter_cdna <- dlply(assemblys, .(accession), function(x){
  print(paste(x$filter1_cdna_tsv))
  tsv = fread(input = as.character(x$filter1_cdna_tsv))
  setnames(tsv, gsub("%","",names(tsv)))
  tsv
})
#
plots_tqcover_cdna_before = list()
for (i in 1:length(dd$accession)) {
  acc = dd$accession[i]
  #
  incontigs_llist = assemblys_contigs_after_cds_llist[[acc]]
  tsv = assemblys_filter_cdna[[acc]]
  no_hits = setdiff(incontigs_llist$contig_name, tsv$qName)
  title = paste(
    "cDNA:",dd$accession[i],dd$tissue[i],"\n",length(no_hits),"out of",
    length(incontigs_llist$contig_name), "contigs were not aligned"
  )
  #
  plots_tqcover_cdna_before[[i]] = ggplot(data = tsv, aes(x = tCoverage, y = qCoverage, 
        colour = identity)) +geom_jitter(alpha = 0.5) + geom_density2d() + theme(legend.position = "bottom") +
    ggtitle(title) + scale_x_continuous(limits=c(0,100)) + scale_y_continuous(limits=c(0,100)) +
    scale_colour_gradientn(
      name = "Identity:  ",limits = c(50,100),
      colours = c("red","yellow","green","lightblue","darkblue"),
      breaks = c(50,75,100),labels = c("low(50)","medium(75)","high(100)"),
      guide = guide_colorbar(
        title.theme = element_text(size = 14, angle = 0),title.vjust = 1,
        barheight = 0.6, barwidth = 6, label.theme = element_text(size = 10, angle = 0)
      )
    )
}
Cairo(
  width = 1200, height = 1200,
  file = "cdna_alignment.pdf", type = "pdf", pointsize = 20,
  bg = "transparent", canvas = "white", units = "px", dpi = 60
)
do.call(grid.arrange,  plots_tqcover_cdna_before)
dev.off()
#
# Get the out filter contigs LLIST
#
assemblys$contigs_after_cdna_llist <-
  daply(assemblys, .(accession), function(x) {
    find_file(x$fullpath, "contigs_after_cdna.llist", recursive = T)
  })
assemblys_contigs_after_cdna_llist <- dlply(assemblys, .(accession), function(x){
  print(paste(x$contigs_after_cdna_llist))
  list = fread(input = as.character(x$contigs_after_cdna_llist), header = F)
  setnames(list, c("contig_name", "length"))
  list
})
#
# plot the out CDS filter figure
#
plots_tqcover_cdna_after = list()
for (i in 1:length(dd$accession)) {
  acc = dd$accession[i]
  #
  incontigs_llist = assemblys_contigs_after_cds_llist[[acc]]
  tsv = assemblys_filter_cdna[[acc]]
  no_hits = setdiff(incontigs_llist$contig_name, tsv$qName)
  outcontigs_llist = assemblys_contigs_after_cdna_llist[[acc]]
  tsv = filter(assemblys_filter_cds[[acc]], qName %in% outcontigs_llist$contig_name)
  title = paste("filtered cDNA:",dd$accession[i],dd$tissue[i],"\n",length(no_hits),"out of",
                length(incontigs_llist$contig_name), "contigs were not aligned")
  #
  plots_tqcover_cdna_after[[i]] = ggplot(data = tsv, aes(x = tCoverage, y = qCoverage, 
       colour = identity)) + geom_jitter(alpha = 0.5) + geom_density2d() + theme(legend.position = "bottom") +
    ggtitle(title) + scale_x_continuous(limits=c(0,100)) + scale_y_continuous(limits=c(0,100)) +
    scale_colour_gradientn(
      name = "Identity:  ",limits = c(50,100),
      colours = c("red","yellow","green","lightblue","darkblue"),
      breaks = c(50,75,100),labels = c("low(50)","medium(75)","high(100)"),
      guide = guide_colorbar(
        title.theme = element_text(size = 14, angle = 0),title.vjust = 1,
        barheight = 0.6, barwidth = 6, label.theme = element_text(size = 10, angle = 0)
      )
    )
}
Cairo(
  width = 1200, height = 1200,
  file = "cdna_alignment_filtered.pdf", type = "pdf", pointsize = 20,
  bg = "transparent", canvas = "white", units = "px", dpi = 60
)
do.call(grid.arrange,  plots_tqcover_cdna_after)
dev.off()
#
##################################################################################
####### REFSeq PEP ###############################################################
##################################################################################
#
assemblys$filter2_refseqpep_tsv <-
  daply(assemblys, .(accession), function(x) {
    find_file(x$fullpath, "blat_cdna_refseq_pep.best.tsv", recursive = T)
  })
assemblys_filter_refseq_pep <- dlply(assemblys, .(accession), function(x){
  print(paste(x$filter2_refseqpep_tsv))
  tsv = fread(input = as.character(x$filter2_refseqpep_tsv))
  setnames(tsv, gsub("%","",names(tsv)))
  tsv
})
#
plots_tqcover_refpep_before = list()
for (i in 1:length(dd$accession)) {
  acc = dd$accession[i]
  #
  incontigs_llist = assemblys_contigs_after_cdna_llist[[acc]]
  tsv = assemblys_filter_refseq_pep[[acc]]
  no_hits = setdiff(incontigs_llist$contig_name, tsv$tName)
  title = paste(
    "refSEQ PEP:",dd$accession[i],dd$tissue[i],"\n",length(no_hits),"out of",
    length(incontigs_llist$contig_name), "contigs were not aligned"
  )
  #
  plots_tqcover_refpep_before[[i]] = ggplot(data = tsv, aes(x = qCoverage, y = tCoverage, 
    colour = identity)) +geom_jitter(alpha = 0.5) + geom_density2d() + theme(legend.position = "bottom") +
    ggtitle(title) + scale_x_continuous(limits=c(0,100)) + scale_y_continuous(limits=c(0,100)) +
    scale_colour_gradientn(
      name = "Identity:  ",limits = c(50,100),
      colours = c("red","yellow","green","lightblue","darkblue"),
      breaks = c(50,75,100),labels = c("low(50)","medium(75)","high(100)"),
      guide = guide_colorbar(
        title.theme = element_text(size = 14, angle = 0),title.vjust = 1,
        barheight = 0.6, barwidth = 6, label.theme = element_text(size = 10, angle = 0)
      )
    )
}
Cairo(
  width = 1200, height = 1200,
  file = "refseq_pep_alignment.pdf", type = "pdf", pointsize = 20,
  bg = "transparent", canvas = "white", units = "px", dpi = 60
)
do.call(grid.arrange,  plots_tqcover_refpep_before)
dev.off()
#
# Get the out filter contigs LLIST
#
assemblys$contigs_after_refseq_pep_llist <-
  daply(assemblys, .(accession), function(x) {
    find_file(x$fullpath, "contigs_after_refseq_pep.llist", recursive = T)
  })
assemblys_contigs_after_refseq_pep_llist <- dlply(assemblys, .(accession), function(x){
  print(paste(x$contigs_after_refseq_pep_llist))
  list = fread(input = as.character(x$contigs_after_refseq_pep_llist), header = F)
  setnames(list, c("contig_name", "length"))
  list
})
#
# plot the out CDS filter figure
#
plots_tqcover_refseq_pep_after = list()
for (i in 1:length(dd$accession)) {
  acc = dd$accession[i]
  #
  incontigs_llist = assemblys_contigs_after_cds_llist[[acc]]
  tsv = assemblys_filter_refseq_pep[[acc]]
  no_hits = setdiff(incontigs_llist$contig_name, tsv$tName)
  outcontigs_llist = assemblys_contigs_after_refseq_pep_llist[[acc]]
  tsv = filter(tsv, tName %in% outcontigs_llist$contig_name)
  title = paste("filtered refSeq PEP:",dd$accession[i],dd$tissue[i],"\n",length(no_hits),"out of",
             length(incontigs_llist$contig_name), "contigs were not aligned")
  #
  plots_tqcover_refseq_pep_after[[i]] = ggplot(data = tsv, aes(x = tCoverage, y = qCoverage, 
    colour = identity)) + geom_jitter(alpha = 0.5) + geom_density2d() + theme(legend.position = "bottom") +
    ggtitle(title) + scale_x_continuous(limits=c(0,100)) + scale_y_continuous(limits=c(0,100)) +
    scale_colour_gradientn(
      name = "Identity:  ",limits = c(50,100),
      colours = c("red","yellow","green","lightblue","darkblue"),
      breaks = c(50,75,100),labels = c("low(50)","medium(75)","high(100)"),
      guide = guide_colorbar(
        title.theme = element_text(size = 14, angle = 0),title.vjust = 1,
        barheight = 0.6, barwidth = 6, label.theme = element_text(size = 10, angle = 0)
      )
    )
}
Cairo(
  width = 1200, height = 1200,
  file = "refseq_pep_alignment_filtered.pdf", type = "pdf", pointsize = 20,
  bg = "transparent", canvas = "white", units = "px", dpi = 60
)
do.call(grid.arrange,  plots_tqcover_refseq_pep_after)
dev.off()