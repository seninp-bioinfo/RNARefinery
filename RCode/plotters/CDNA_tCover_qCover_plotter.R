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
percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}
#
# set folders
#
archive_dir <-
  "/work2/project/sigenae/Project_RNA_Refinery/assembled"
setwd(archive_dir)
#
# FUNCTION: list all the folders within the working directory
#
list_dirs <-
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
dirs <- list_dirs(archive_dir, full.names = T)
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
df = ldply(assemblys_contigs_llist,function(x){length(x$contig_name)})
setnames(df,c("accession","contigs"))
assemblys=merge(assemblys,df)
#
#
arranged = dlply(assemblys, .(tissue), function(x){x})
plots_assemblys_by_tissue = llply(arranged, function(x){
  tissue = unique(x$tissue)
  plot_assembly_results <- ggplot(data=x, aes(x=reads, y=contigs, color=accession)) + 
    theme_bw() + geom_jitter(lwd=2, alpha=0.7) + guides(color=guide_legend(nrow=4,byrow=TRUE)) + 
    ggtitle(tissue) + theme(legend.position="bottom", legend.text=element_text(size=8),
        legend.title=element_blank(), legend.key.size = unit(0.3, "lines"))
  plot_assembly_results
})
Cairo(width = 1300, height = 1300, file = "rnarefinery_plots_assemblys_by_tissue.pdf", type = "pdf", pointsize = 20,
      bg = "transparent", canvas = "white", units = "px", dpi = 60 )
do.call(grid.arrange,  plots_assemblys_by_tissue)
dev.off()
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

# get the TSV filenames
#
tsv_mask <- "blat_raw_cds.best.tsv"
assemblys$filter0_cds_tsv <- daply(assemblys, .(accession), function(x) {
    find_file(x$fullpath, tsv_mask, recursive = T)  })

# get the TSV data
#
cname <- "filter0_cds_tsv"
assemblys_filter_cds_tsv <- dlply(assemblys, .(accession), function(x){
  print(paste(x[,paste(cname)]))
  tsv = fread(input = as.character(x[,paste(cname)]))
  setnames(tsv, gsub("%","",names(tsv)))
  tsv
})

# plot the CDS alignment
#
plots_tqcover_cds_before = list()
for (i in 1:length(dd$accession)) {
  #
  acc = dd$accession[i]
  #
  incontigs_llist = assemblys_contigs_llist[[acc]]
  tsv = assemblys_filter_cds_tsv[[acc]]
  #
  no_hits = setdiff(incontigs_llist$contig_name, tsv$qName)
  title = paste("CDS:",dd$accession[i],", ", dd$tissue[i],"\n",length(no_hits),
    " out of ",length(incontigs_llist$contig_name), " contigs were not aligned (", 
    percent(length(no_hits)/length(incontigs_llist$contig_name)),")",sep="")
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

# Get the LLIST filenames
#
llist_mask <- "contigs_after_cds.llist"
assemblys$contigs_after_cds_llist <-
  daply(assemblys, .(accession), function(x) {
    find_file(x$fullpath, llist_mask, recursive = T)
  })

# get the LLIST data
#
cname <- "contigs_after_cds_llist"
assemblys_contigs_after_cds_llist <- dlply(assemblys, .(accession), function(x){
  print(paste(x[,paste(cname)]))
  list = fread(input = as.character(x[,paste(cname)]), header = F)
  setnames(list, c("contig_name", "length"))
  list
})

# plot the out CDS filter figure
#
plots_tqcover_cds_after = list()
for (i in 1:length(dd$accession)) {
  acc = dd$accession[i]
  #
  incontigs_llist = assemblys_contigs_llist[[acc]]
  tsv = assemblys_filter_cds_tsv[[acc]]
  outcontigs_llist = assemblys_contigs_after_cds_llist[[acc]]
  #
  no_hits = setdiff(incontigs_llist$contig_name, tsv$qName)
  tsv_out = filter(tsv, qName %in% outcontigs_llist$contig_name)
  title = paste("filtered CDS: ",dd$accession[i], ", ", dd$tissue[i],"\nin ",length(incontigs_llist$contig_name),
    ", out ",length(outcontigs_llist$contig_name)," (no hits: ",length(no_hits), ", excluded: ",
     length(tsv$qName) - (length(outcontigs_llist$contig_name) - length(no_hits)), ")", sep="")
  #
  plots_tqcover_cds_after[[i]] = ggplot(data = tsv_out, aes(x = tCoverage, y = qCoverage, 
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

# get the TSV filenames
#
tsv_mask <- "blat_cds_cdna.best.tsv"
assemblys$filter1_cdna_tsv <- daply(assemblys, .(accession), function(x) {
  find_file(x$fullpath, tsv_mask, recursive = T)  })

# get the TSV data
#
cname <- "filter1_cdna_tsv"
assemblys_filter_cdna_tsv <- dlply(assemblys, .(accession), function(x){
  print(paste(x[,paste(cname)]))
  tsv = fread(input = as.character(x[,paste(cname)]))
  setnames(tsv, gsub("%","",names(tsv)))
  tsv
})

# plot the CDNA alignment
#
plots_tqcover_cdna_before = list()
for (i in 1:length(dd$accession)) {
  #
  acc = dd$accession[i]
  #
  incontigs_llist = assemblys_contigs_after_cds_llist[[acc]]
  tsv = assemblys_filter_cdna_tsv[[acc]]
  #
  no_hits = setdiff(incontigs_llist$contig_name, tsv$qName)
  title = paste("CDNA:",dd$accession[i],", ",dd$tissue[i],"\n",length(no_hits)," out of ",
                length(incontigs_llist$contig_name), " contigs were not aligned (", 
                percent(length(no_hits)/length(incontigs_llist$contig_name)),")",sep="")
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

# Get the LLIST filenames
#
llist_mask <- "contigs_after_cdna.llist"
assemblys$contigs_after_cdna_llist <-
  daply(assemblys, .(accession), function(x) {
    find_file(x$fullpath, llist_mask, recursive = T)
  })

# get the LLIST data
#
cname <- "contigs_after_cdna_llist"
assemblys_contigs_after_cdna_llist <- dlply(assemblys, .(accession), function(x){
  print(paste(x[,paste(cname)]))
  list = fread(input = as.character(x[,paste(cname)]), header = F)
  setnames(list, c("contig_name", "length"))
  list
})

# plot the out CDNA filter figure
#
plots_tqcover_cdna_after = list()
for (i in 1:length(dd$accession)) {
  acc = dd$accession[i]
  #
  incontigs_llist = assemblys_contigs_after_cds_llist[[acc]]
  tsv = assemblys_filter_cdna_tsv[[acc]]
  outcontigs_llist = assemblys_contigs_after_cdna_llist[[acc]]
  #
  no_hits = setdiff(incontigs_llist$contig_name, tsv$qName)
  tsv_out = filter(tsv, qName %in% outcontigs_llist$contig_name)
  title = paste("filtered CDNA: ",dd$accession[i], ", ", dd$tissue[i],"\nin ",length(incontigs_llist$contig_name),
    ", out ",length(outcontigs_llist$contig_name)," (no hits: ",length(no_hits), ", excluded: ",
     length(tsv$qName) - (length(outcontigs_llist$contig_name) - length(no_hits)), ")", sep="")
  #
  plots_tqcover_cdna_after[[i]] = ggplot(data = tsv_out, aes(x = tCoverage, y = qCoverage, 
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
Cairo(width = 1200, height = 1200, file = "cdna_alignment_filtered.pdf", type = "pdf", pointsize = 20,
  bg = "transparent", canvas = "white", units = "px", dpi = 60 )
do.call(grid.arrange,  plots_tqcover_cdna_after)
dev.off()
#
##################################################################################
####### REFSeq PEP ###############################################################
##################################################################################

# get the TSV filenames
#
tsv_mask <- "blat_cdna_refseq_pep.best.tsv"
assemblys$filter2_refseqpep_tsv <- daply(assemblys, .(accession), function(x) {
  find_file(x$fullpath, tsv_mask, recursive = T)  })

# get the TSV data
#
cname <- "filter2_refseqpep_tsv"
assemblys_filter_refseq_pep_tsv <- dlply(assemblys, .(accession), function(x){
  print(paste(x[,paste(cname)]))
  tsv = fread(input = as.character(x[,paste(cname)]))
  setnames(tsv, gsub("%","",names(tsv)))
  tsv
})

# plot the REFSEQ PEP alignment
#
plots_tqcover_pep_before = list()
for (i in 1:length(dd$accession)) {
  #
  acc = dd$accession[i]
  #
  incontigs_llist = assemblys_contigs_after_cdna_llist[[acc]]
  tsv = assemblys_filter_refseq_pep_tsv[[acc]]
  #
  no_hits = setdiff(incontigs_llist$contig_name, tsv$tName)
  title = paste("REFSeq PEP:",dd$accession[i],", ",dd$tissue[i],"\n",length(no_hits)," out of ",
    length(incontigs_llist$contig_name), " contigs were not aligned (", 
    percent(length(no_hits)/length(incontigs_llist$contig_name)),")",sep="")
  #
  plots_tqcover_pep_before[[i]] = ggplot(data = tsv, aes(x = qCoverage, y = tCoverage, 
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
do.call(grid.arrange,  plots_tqcover_pep_before)
dev.off()

# Get the LLIST filenames
#
llist_mask <- "contigs_after_refseq_pep.llist"
assemblys$contigs_after_refseq_pep_llist <-
  daply(assemblys, .(accession), function(x) {
    find_file(x$fullpath, llist_mask, recursive = T)
  })

# get the LLIST data
#
cname <- "contigs_after_refseq_pep_llist"
assemblys_contigs_after_refseq_pep_llist <- dlply(assemblys, .(accession), function(x){
  print(paste(x[,paste(cname)]))
  list = fread(input = as.character(x[,paste(cname)]), header = F)
  setnames(list, c("contig_name", "length"))
  list
})

# plot the out REFSEQ PEP filter figure
#
plots_tqcover_refseq_pep_after = list()
for (i in 1:length(dd$accession)) {
  acc = dd$accession[i]
  #
  incontigs_llist = assemblys_contigs_after_cdna_llist[[acc]]
  tsv = assemblys_filter_refseq_pep_tsv[[acc]]
  outcontigs_llist = assemblys_contigs_after_refseq_pep_llist[[acc]]
  #
  no_hits = setdiff(incontigs_llist$contig_name, tsv$tName)
  tsv_out = filter(tsv, tName %in% outcontigs_llist$contig_name)
  title = paste("filtered REFSeq PEP: ",dd$accession[i], ", ", dd$tissue[i],"\nin ",length(incontigs_llist$contig_name),
    ", out ",length(outcontigs_llist$contig_name)," (no hits: ",length(no_hits), ", excluded: ",
     length(tsv$tName) - (length(outcontigs_llist$contig_name) - length(no_hits)), ")", sep="")
  #
  plots_tqcover_refseq_pep_after[[i]] = ggplot(data = tsv_out, aes(x = qCoverage, y = tCoverage, 
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
Cairo(width = 1200, height = 1200, file = "refseq_pep_alignment_filtered.pdf", type = "pdf", pointsize = 20,
  bg = "transparent", canvas = "white", units = "px", dpi = 60 )
do.call(grid.arrange,  plots_tqcover_refseq_pep_after)
dev.off()
##################################################################################
####### REFSeq DNA ###############################################################
##################################################################################

# get the TSV filenames
#
tsv_mask <- "blat_pep_refseq_dna.best.tsv"
assemblys$filter3_refseq_dna_tsv <- daply(assemblys, .(accession), function(x) {
  find_file(x$fullpath, tsv_mask, recursive = T)  })

# get the TSV data
#
cname <- "filter3_refseq_dna_tsv"
assemblys_filter_refseq_dna_tsv <- dlply(assemblys, .(accession), function(x){
  print(paste(x[,paste(cname)]))
  tsv = fread(input = as.character(x[,paste(cname)]))
  setnames(tsv, gsub("%","",names(tsv)))
  tsv
})

# plot the CDNA alignment
#
plots_tqcover_refseq_dna_before = list()
for (i in 1:length(dd$accession)) {
  #
  acc = dd$accession[i]
  #
  incontigs_llist = assemblys_contigs_after_refseq_pep_llist[[acc]]
  tsv = assemblys_filter_refseq_dna_tsv[[acc]]
  #
  no_hits = setdiff(incontigs_llist$contig_name, tsv$qName)
  title = paste("REFSeq DNA:",dd$accession[i],", ",dd$tissue[i],"\n",length(no_hits)," out of ",
                length(incontigs_llist$contig_name), " contigs were not aligned (", 
    percent(length(no_hits)/length(incontigs_llist$contig_name)),")",sep="")
  #
  plots_tqcover_refseq_dna_before[[i]] = ggplot(data = tsv, aes(x = tCoverage, y = qCoverage, 
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
  file = "refseq_dna_alignment.pdf", type = "pdf", pointsize = 20,
  bg = "transparent", canvas = "white", units = "px", dpi = 60
)
do.call(grid.arrange,  plots_tqcover_refseq_dna_before)
dev.off()

# Get the LLIST filenames
#
llist_mask <- "contigs_after_refseq_dna.llist"
assemblys$contigs_after_refseq_dna_llist <-
  daply(assemblys, .(accession), function(x) {
    find_file(x$fullpath, llist_mask, recursive = T)
  })

# get the LLIST data
#
cname <- "contigs_after_refseq_dna_llist"
assemblys_contigs_after_refseq_dna_llist <- dlply(assemblys, .(accession), function(x){
  print(paste(x[,paste(cname)]))
  list = fread(input = as.character(x[,paste(cname)]), header = F)
  setnames(list, c("contig_name", "length"))
  list
})

# plot the out CDNA filter figure
#
plots_tqcover_refseq_dna_after = list()
for (i in 1:length(dd$accession)) {
  acc = dd$accession[i]
  #
  incontigs_llist = assemblys_contigs_after_refseq_pep_llist[[acc]]
  tsv = assemblys_filter_refseq_dna_tsv[[acc]]
  outcontigs_llist = assemblys_contigs_after_refseq_dna_llist[[acc]]
  #
  no_hits = setdiff(incontigs_llist$contig_name, tsv$qName)
  tsv_out = filter(tsv, qName %in% outcontigs_llist$contig_name)
  title = paste("filtered REFSeq DNA: ",dd$accession[i], ", ", dd$tissue[i],"\nin ",length(incontigs_llist$contig_name),
    ", out ",length(outcontigs_llist$contig_name)," (no hits: ",length(no_hits), ", excluded: ",
     length(tsv$qName) - (length(outcontigs_llist$contig_name) - length(no_hits)), ")", sep="")
  #
  plots_tqcover_refseq_dna_after[[i]] = ggplot(data = tsv_out, aes(x = tCoverage, y = qCoverage, 
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
Cairo(width = 1200, height = 1200, file = "refseq_dna_alignment_filtered.pdf", type = "pdf", pointsize = 20,
  bg = "transparent", canvas = "white", units = "px", dpi = 60 )
do.call(grid.arrange,  plots_tqcover_refseq_dna_after)
dev.off()
##################################################################################
####### GENOME ###################################################################
##################################################################################

# get the TSV filenames
#
tsv_mask <- "blat_pep_genome.best.tsv"
assemblys$filter3_genome_tsv <- daply(assemblys, .(accession), function(x) {
  find_file(x$fullpath, tsv_mask, recursive = T)  })

# get the TSV data
#
cname <- "filter3_genome_tsv"
assemblys_filter_genome_tsv <- dlply(assemblys, .(accession), function(x){
  print(paste(x[,paste(cname)]))
  tsv = fread(input = as.character(x[,paste(cname)]))
  setnames(tsv, gsub("%","",names(tsv)))
  tsv
})

# plot the GENOME alignment
#
plots_tqcover_genome_before = list()
for (i in 1:length(dd$accession)) {
  #
  acc = dd$accession[i]
  #
  incontigs_llist = assemblys_contigs_after_refseq_pep_llist[[acc]]
  tsv = assemblys_filter_genome_tsv[[acc]]
  #
  no_hits = setdiff(incontigs_llist$contig_name, tsv$qName)
  title = paste("GENOME:",dd$accession[i],", ",dd$tissue[i],"\n",length(no_hits)," out of ",
                length(incontigs_llist$contig_name), " contigs were not aligned",sep="")
  #
  plots_tqcover_genome_before[[i]] = ggplot(data = tsv, aes(x = tCoverage, y = qCoverage, 
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
  file = "refseq_pep_genome_alignment.pdf", type = "pdf", pointsize = 20,
  bg = "transparent", canvas = "white", units = "px", dpi = 60
)
do.call(grid.arrange,  plots_tqcover_genome_before)
dev.off()

# Get the LLIST filenames
#
llist_mask <- "contigs_after_genome.llist"
assemblys$contigs_after_genome_llist <-
  daply(assemblys, .(accession), function(x) {
    find_file(x$fullpath, llist_mask, recursive = T)
  })

# get the LLIST data
#
cname <- "contigs_after_genome_llist"
assemblys_contigs_after_genome_llist <- dlply(assemblys, .(accession), function(x){
  print(paste(x[,paste(cname)]))
  list = fread(input = as.character(x[,paste(cname)]), header = F)
  setnames(list, c("contig_name", "length"))
  list
})

# plot the out GENOME filter figure
#
plots_tqcover_genome_after = list()
for (i in 1:length(dd$accession)) {
  acc = dd$accession[i]
  #
  incontigs_llist = assemblys_contigs_after_refseq_pep_llist[[acc]]
  tsv = assemblys_filter_genome_tsv[[acc]]
  outcontigs_llist = assemblys_contigs_after_genome_llist[[acc]]
  #
  no_hits = setdiff(incontigs_llist$contig_name, tsv$qName)
  tsv_out = filter(tsv, qName %in% outcontigs_llist$contig_name)
  title = paste("filtered GENOME: ",dd$accession[i], ", ", dd$tissue[i],"\nin ",length(incontigs_llist$contig_name),
    ", out ",length(outcontigs_llist$contig_name)," (no hits: ",length(no_hits), ", excluded: ",
     length(tsv$qName) - (length(outcontigs_llist$contig_name) - length(no_hits)), ")", sep="")
  #
  plots_tqcover_genome_after[[i]] = ggplot(data = tsv_out, aes(x = tCoverage, y = qCoverage, 
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
Cairo(width = 1200, height = 1200, file = "refseq_pep_genome_alignment_filtered.pdf", type = "pdf", pointsize = 20,
  bg = "transparent", canvas = "white", units = "px", dpi = 60 )
do.call(grid.arrange,  plots_tqcover_genome_after)
dev.off()

##################################################################################
####### GENOME AFTER REFSEQ DNA ##################################################
##################################################################################

# get the TSV filenames
#
tsv_mask <- "blat_reseq_dna_genome.best.tsv"
assemblys$filter4_genome_tsv <- daply(assemblys, .(accession), function(x) {
  find_file(x$fullpath, tsv_mask, recursive = T)  })

# get the TSV data
#
cname <- "filter4_genome_tsv"
assemblys_filter_genome2_tsv <- dlply(assemblys, .(accession), function(x){
  print(paste(x[,paste(cname)]))
  tsv = fread(input = as.character(x[,paste(cname)]))
  setnames(tsv, gsub("%","",names(tsv)))
  tsv
})

# plot the GENOME alignment
#
plots_tqcover_genome_before = list()
for (i in 1:length(dd$accession)) {
  #
  acc = dd$accession[i]
  #
  incontigs_llist = assemblys_contigs_after_refseq_dna_llist[[acc]]
  tsv = assemblys_filter_genome2_tsv[[acc]]
  #
  no_hits = setdiff(incontigs_llist$contig_name, tsv$qName)
  title = paste("GENOME2:",dd$accession[i],", ",dd$tissue[i],"\n",length(no_hits)," out of ",
                length(incontigs_llist$contig_name), " contigs were not aligned",sep="")
  #
  plots_tqcover_genome_before[[i]] = ggplot(data = tsv, aes(x = tCoverage, y = qCoverage, 
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
  file = "refseq_dna_genome_alignment.pdf", type = "pdf", pointsize = 20,
  bg = "transparent", canvas = "white", units = "px", dpi = 60
)
do.call(grid.arrange,  plots_tqcover_genome_before)
dev.off()

# Get the LLIST filenames
#
llist_mask <- "contigs_after_genome2.llist"
assemblys$contigs_after_genome2_llist <-
  daply(assemblys, .(accession), function(x) {
    find_file(x$fullpath, llist_mask, recursive = T)
  })

# get the LLIST data
#
cname <- "contigs_after_genome2_llist"
assemblys_contigs_after_genome2_llist <- dlply(assemblys, .(accession), function(x){
  print(paste(x[,paste(cname)]))
  list = fread(input = as.character(x[,paste(cname)]), header = F)
  setnames(list, c("contig_name", "length"))
  list
})

# plot the out GENOME filter figure
#
plots_tqcover_genome_after = list()
for (i in 1:length(dd$accession)) {
  acc = dd$accession[i]
  #
  incontigs_llist = assemblys_contigs_after_refseq_dna_llist[[acc]]
  tsv = assemblys_filter_genome2_tsv[[acc]]
  outcontigs_llist = assemblys_contigs_after_genome2_llist[[acc]]
  #
  no_hits = setdiff(incontigs_llist$contig_name, tsv$qName)
  tsv_out = filter(tsv, qName %in% outcontigs_llist$contig_name)
  title = paste("filtered GENOME2: ",dd$accession[i], ", ", dd$tissue[i],"\nin ",length(incontigs_llist$contig_name),
    ", out ",length(outcontigs_llist$contig_name)," (no hits: ",length(no_hits), ", excluded: ",
     length(tsv$qName) - (length(outcontigs_llist$contig_name) - length(no_hits)), ")", sep="")
  #
  plots_tqcover_genome_after[[i]] = ggplot(data = tsv_out, aes(x = tCoverage, y = qCoverage, 
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
Cairo(width = 1200, height = 1200, file = "refseq_dna_genome_alignment_filtered.pdf", type = "pdf", pointsize = 20,
  bg = "transparent", canvas = "white", units = "px", dpi = 60 )
do.call(grid.arrange,  plots_tqcover_genome_after)
dev.off()
#
#
#
dd=data.frame(accession=names(assemblys_contigs_llist), stringsAsFactors=F)
dd$raw = daply(dd,.(accession),function(x){length(assemblys_contigs_llist[[x$accession]]$contig_name)})
dd$cds = daply(dd,.(accession),function(x){length(assemblys_contigs_after_cds_llist[[x$accession]]$contig_name)})
dd$cdna = daply(dd,.(accession),function(x){length(assemblys_contigs_after_cdna_llist[[x$accession]]$contig_name)})
dd$refseq_pep = daply(dd,.(accession),function(x){length(assemblys_contigs_after_refseq_pep_llist[[x$accession]]$contig_name)})
dd$genome = daply(dd,.(accession),function(x){length(assemblys_contigs_after_genome_llist[[x$accession]]$contig_name)})
dd$refseq_dna = daply(dd,.(accession),function(x){length(assemblys_contigs_after_refseq_dna_llist[[x$accession]]$contig_name)})
dd$genome2 = daply(dd,.(accession),function(x){length(assemblys_contigs_after_genome2_llist[[x$accession]]$contig_name)})
dd$tissue = daply(dd,.(accession),function(x){assemblys[assemblys$accession==x$accession,]$tissue})

df = melt(dd[,1:6],id.vars=c("accession"))
p1 = ggplot(df, aes(x=variable,y=value,group=accession)) + geom_line() + ggtitle("Filtering option 1")

df2 = melt(dd[,c(1,2,3,4,5,7,8)],id.vars=c("accession"))
p2 = ggplot(df2, aes(x=variable,y=value,group=accession)) + geom_line() + ggtitle("Filtering option 2 (via REFSeq DNA)")
Cairo(width = 900, height = 600, file = "filtering_process.pdf", type = "pdf", pointsize = 20,
      bg = "transparent", canvas = "white", units = "px", dpi = 60 )
grid.arrange(p1, p2, ncol=2)
dev.off()

arranged = dlply(dd, .(tissue), function(x){x})
plots_filter_by_tissue = llply(arranged, function(x){
  tissue = unique(x$tissue)
  df = melt(x[,c(1,2,3,4,5,7,8)],id.vars=c("accession"))
  plot_assembly_results <- ggplot(df, aes(x=variable,y=value,group=accession)) + geom_line() + 
    ggtitle(tissue) + theme(legend.position="bottom", legend.text=element_text(size=8),
                            legend.title=element_blank(), legend.key.size = unit(0.3, "lines"))
  plot_assembly_results
})

Cairo(width = 1300, height = 1300, file = "filter_by_tissue2.pdf", type = "pdf", pointsize = 20,
        bg = "transparent", canvas = "white", units = "px", dpi = 60 )
do.call(grid.arrange,  plots_filter_by_tissue)
dev.off()