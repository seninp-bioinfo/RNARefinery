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
""
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
  find_complete(paste(x$fullpath,"/filters/",sep=""),"contigs_after_cdna.llist")
})
#
#
# attach the run tissue description
#
properties = read.table("/home/psenin/.rnarefinery/db.properties", header = F, as.is = T)
db_credentials = data.frame(t(properties$V2), stringsAsFactors = F); names(db_credentials) = properties$V1
session <-
  dbConnect(
    MySQL(), host = db_credentials$host, dbname = db_credentials$db,
    user = db_credentials$user, password = db_credentials$password
  )
a_tissue <- ddply(assemblys, .(accession), function(x) {
  tissue = dbGetQuery(session, paste(
    "select tissue from sra_runs where run_accession=", shQuote(x$accession)
  ))
  tissue
})
names(a_tissue) <- c("accession","tissue")
assemblys = merge(assemblys,a_tissue)
#
#
# filter out only complete assemblies
assemblys <- assemblys[grep("^complete$", assemblys$complete),]
#
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
# get the first CDS alignment
#
assemblys$filter0_cds_tsv <-
  daply(assemblys, .(accession), function(x) {
    find_file(x$fullpath, "blat_raw_cds.best.tsv", recursive = T)
  })

assemblys_filter_cds <- dlply(assemblys, .(accession), function(x){
  print(paste(x$filter0_cds_tsv))
  tsv = fread(input = as.character(x$filter0_cds_tsv))
  tsv
})
#
# get the second CDNA alignment
#
assemblys$filter1_cdna_tsv <-
  daply(assemblys, .(accession), function(x) {
    find_file(x$fullpath, "blat_cds_cdna.best.tsv", recursive = T)
  })
assemblys_filter_cdna <- dlply(assemblys, .(accession), function(x){
  print(paste(x$filter1_cdna_tsv))
  tsv = fread(input = as.character(x$filter1_cdna_tsv))
  tsv
})
#
# get the refseq PEP alignment
#
assemblys$filter2_refseqpep_tsv <-
  daply(assemblys, .(accession), function(x) {
    find_file(x$fullpath, "blat_cdna_refseq_pep.best.tsv", recursive = T)
  })
assemblys_filter_cdna <- dlply(assemblys, .(accession), function(x){
  print(paste(x$filter2_refseqpep_tsv))
  tsv = fread(input = as.character(x$filter2_refseqpep_tsv))
  tsv
})

# get the contigs after 75 squared
#
assemblys$squared75_llist <-
  daply(assemblys, .(accession), function(x) {
    find_file(x$fullpath, "tcover_contigs_after_cds.llist", recursive = T)
  })
squared_contigs_llist <- dlply(assemblys, .(accession), function(x){
  print(paste(x$squared75_llist))
  list = fread(input = as.character(x$squared75_llist), header = F)
  setnames(list, c("contig_name", "length"))
  list
})

plots_tqcover = list()
for(i in 1:length(dd$accession)){
  acc = dd$accession[i]
  dataset = assemblys_filter_cds[[acc]]
  outcontigs_llist = squared_contigs_llist[[acc]]
  setnames(dataset, gsub("%","",names(dataset)))
  #dataset = filter(dataset,qLength>500)
  dataset = dataset[dataset$qName %in% outcontigs_llist$contig_name,]
  plots_tqcover[[i]] = ggplot(data = dataset, aes(x = tCoverage, y = qCoverage, colour=identity)) + 
    geom_jitter(alpha=0.5) + geom_density2d() + theme(legend.position="bottom") +
    ggtitle(paste("tqCoverage (q>500)",paste(dd$accession[i],dd$tissue[i]))) +
    scale_colour_gradientn(name = "Identity:  ",limits=c(50,100),
    colours=c("red","yellow","green","lightblue","darkblue"),
    breaks=c(50,75,100),labels=c("low","medium","high"),
    guide = guide_colorbar(title.theme=element_text(size=14, angle=0),title.vjust=1,
    barheight=0.6, barwidth=6, label.theme=element_text(size=10, angle=0)))
}
Cairo(width = 1600, height = 1600, 
      file="after_75cutoff_tcoverage.pdf", type="pdf", pointsize=24, 
      bg = "transparent", canvas = "white", units = "px", dpi = 55)
do.call(grid.arrange,  plots_tqcover)
dev.off()




# OVARY and others
#
#
#
idx <- c(grep("ovary",assemblys$tissue), grep("magnum",assemblys$tissue)[1:2], 
         grep("testis",assemblys$tissue)[1:2])
#
# pull these out
#
dd = assemblys[idx,]
dd = arrange(dd,tissue)
#

#
# plot assemblys on SPOTS - CONTIGS axis
#
which(names(assemblys_contigs_llist) %in% dd$accession)
df = ldply(assemblys_contigs_llist,function(x){length(x$contig_name)})[df$accession %in% dd$accession,]
df$tissue = daply(df, .(accession), function(x){dd[dd$accession==x$accession,]$tissue})
df$reads = daply(df, .(accession), function(x){dd[dd$accession==x$accession,]$reads})
setnames(df,c("accession","contigs","tissue","reads"))
plot_spots_assemblys <- ggplot(data=df,aes(x=reads/1e+6,y=contigs,color=paste(accession,tissue))) + theme_bw() +
  geom_point(size=10) + ggtitle("Sequenced reads (millions) versus contigs")
#
# plot contigs length histograms
#
plots_contigs_len = list()
for(i in 1:length(dd$accession)){
  acc = dd$accession[i]
  dataset = assemblys_contigs_llist[[acc]]
  plots_contigs_len[[i]] = ggplot(data = dataset, aes(x = length)) +
    scale_x_continuous(limits = c(0,5000)) + scale_y_log10() +
    geom_histogram(binwidth = 100) + ggtitle(paste("Contigs len for",paste(dd$accession[i],dd$tissue[i])))
}
#
# plot the QUERY and TEMPLATE overlap histograms
#
plots_qcover_len = list()
plots_tcover_len = list()
plots_tcover_identity = list()
for(i in 1:length(dd$accession)){
  acc = dd$accession[i]
  dataset = assemblys_filter_cds[[acc]]
  setnames(dataset, gsub("%","",names(dataset)))
  plots_qcover_len[[i]] = ggplot(data=dataset, aes(x=qCoverage)) + 
    geom_histogram(binwidth = 5) + ggtitle(paste("qCoverage for",paste(dd$accession[i],dd$tissue[i])))
  plots_tcover_len[[i]] = ggplot(data=dataset, aes(x=tCoverage)) + 
    geom_histogram(binwidth = 5) + ggtitle(paste("tCoverage for",paste(dd$accession[i],dd$tissue[i])))
  plots_tcover_identity[[i]]= ggplot(data=dataset, aes(x=tCoverage, y=identity)) + 
    geom_jitter(alpha=0.6) + geom_density2d() + 
    ggtitle(paste("tCoverage VS identity for",paste(dd$accession[i],dd$tissue[i])))
}
#
#
#
Cairo(width = 2000, height = 1300, 
      file="test2.pdf", type="pdf", pointsize=24, 
      bg = "transparent", canvas = "white", units = "px", dpi = 60)
print(grid.arrange(plot_spots_assemblys,
                   arrangeGrob(plots_contigs_len[[1]], plots_contigs_len[[2]],plots_contigs_len[[3]], 
                               plots_contigs_len[[4]], plots_contigs_len[[5]],plots_contigs_len[[6]], 
                               plots_contigs_len[[7]], ncol=7),
                   arrangeGrob(plots_qcover_len[[1]], plots_qcover_len[[2]],plots_qcover_len[[3]], 
                               plots_qcover_len[[4]], plots_qcover_len[[5]],plots_qcover_len[[6]], 
                               plots_qcover_len[[7]], ncol=7),     
                   arrangeGrob(plots_tcover_len[[1]], plots_tcover_len[[2]],plots_tcover_len[[3]], 
                               plots_tcover_len[[4]], plots_tcover_len[[5]],plots_tcover_len[[6]], 
                               plots_tcover_len[[7]], ncol=7),
                   arrangeGrob(plots_tcover_identity[[1]], plots_tcover_identity[[2]],plots_tcover_identity[[3]], 
                               plots_tcover_identity[[4]], plots_tcover_identity[[5]],plots_tcover_identity[[6]], 
                               plots_tcover_identity[[7]], ncol=7),     
                   nrow=5))
dev.off()


Cairo(width = 1000, height = 1000, 
      file="test3.pdf", type="pdf", pointsize=24, 
      bg = "transparent", canvas = "white", units = "px", dpi = 60)
print(grid.arrange(plot_spots_assemblys,
                   arrangeGrob(plots_contigs_len[[1]], plots_contigs_len[[2]],plots_contigs_len[[3]], 
                               plots_contigs_len[[4]], plots_contigs_len[[5]],plots_contigs_len[[6]], 
                               plots_contigs_len[[7]], ncol=7),
                   arrangeGrob(plots_qcover_len[[1]], plots_qcover_len[[2]],plots_qcover_len[[3]], 
                               plots_qcover_len[[4]], plots_qcover_len[[5]],plots_qcover_len[[6]], 
                               plots_qcover_len[[7]], ncol=7),     
                   arrangeGrob(plots_tcover_len[[1]], plots_tcover_len[[2]],plots_tcover_len[[3]], 
                               plots_tcover_len[[4]], plots_tcover_len[[5]],plots_tcover_len[[6]], 
                               plots_tcover_len[[7]], ncol=7),
                   arrangeGrob(plots_tcover_identity[[1]], plots_tcover_identity[[2]],plots_tcover_identity[[3]], 
                               plots_tcover_identity[[4]], plots_tcover_identity[[5]],plots_tcover_identity[[6]], 
                               plots_tcover_identity[[7]], ncol=7),     
                   nrow=5))
dev.off()


idx <- sample(1:301,30)
dd = assemblys[idx,]
dd = arrange(dd,tissue)
plots_contigs_len = list()
for(i in 1:length(dd$accession)){
  acc = dd$accession[i]
  dataset = assemblys_contigs_llist[[acc]]
  plots_contigs_len[[i]] = ggplot(data = dataset, aes(x = length)) +
    scale_x_continuous(limits = c(0,5000)) + scale_y_log10() +
    geom_histogram(binwidth = 100) + ggtitle(paste("Contigs len for",paste(dd$accession[i],dd$tissue[i])))
}

plots_tqcover = list()
for(i in 1:length(dd$accession)){
  acc = dd$accession[i]
  dataset = assemblys_filter_cds[[acc]]
  setnames(dataset, gsub("%","",names(dataset)))
  dataset = filter(dataset,qLength>500)
  plots_tqcover[[i]] = ggplot(data = dataset, aes(x = tCoverage, y = qCoverage, colour=identity)) + 
    geom_jitter(alpha=0.5) + geom_density2d() + theme(legend.position="bottom") +
    ggtitle(paste("tqCoverage (q>500)",paste(dd$accession[i],dd$tissue[i]))) +
    scale_colour_gradientn(name = "Identity:  ",limits=c(60,100),
    colours=c("red","yellow","green","lightblue","darkblue"),
    breaks=c(60,80,100),labels=c("low","medium","high"),
    guide = guide_colorbar(title.theme=element_text(size=14, angle=0),title.vjust=1,
    barheight=0.6, barwidth=6, label.theme=element_text(size=10, angle=0)))
}
Cairo(width = 1600, height = 1600, 
      file="before_75cutoff_tcoverage.pdf", type="pdf", pointsize=24, 
      bg = "transparent", canvas = "white", units = "px", dpi = 55)
do.call(grid.arrange,  plots_tqcover)
dev.off()

plots_tqcover = list()
for(i in 1:length(dd$accession)){
  acc = dd$accession[i]
  dataset = assemblys_filter_cds[[acc]]
  outcontigs_llist = squared_contigs_llist[[acc]]
  setnames(dataset, gsub("%","",names(dataset)))
  #dataset = filter(dataset,qLength>500)
  dataset = dataset[dataset$qName %in% outcontigs_llist$contig_name,]
  plots_tqcover[[i]] = ggplot(data = dataset, aes(x = tCoverage, y = qCoverage, colour=identity)) + 
    geom_jitter(alpha=0.5) + geom_density2d() + theme(legend.position="bottom") +
    ggtitle(paste("tqCoverage (q>500)",paste(dd$accession[i],dd$tissue[i]))) +
    scale_colour_gradientn(name = "Identity:  ",limits=c(20,100),
    colours=c("red","yellow","green","lightblue","darkblue"),
    breaks=c(60,80,100),labels=c("low","medium","high"),
    guide = guide_colorbar(title.theme=element_text(size=14, angle=0),title.vjust=1,
            barheight=0.6, barwidth=6, label.theme=element_text(size=10, angle=0)))
}
Cairo(width = 1600, height = 1600, 
      file="after_75cutoff_tcoverage.pdf", type="pdf", pointsize=24, 
      bg = "transparent", canvas = "white", units = "px", dpi = 55)
do.call(grid.arrange,  plots_tqcover)
dev.off()




plots_tqcover_cdna = list()
for(i in 1:length(dd$accession)){
  acc = dd$accession[i]
  dataset = assemblys_filter_cdna[[acc]]
  setnames(dataset, gsub("%","",names(dataset)))
  dataset = filter(dataset,qLength>500)
  plots_tqcover_cdna[[i]] = ggplot(data = dataset, aes(x = tCoverage, y = qCoverage, colour=identity)) + 
    geom_jitter(alpha=0.5) + geom_density2d() + theme(legend.position="bottom") +
    ggtitle(paste("tqCoverage (q>500)",paste(dd$accession[i],dd$tissue[i]))) +
    scale_colour_gradientn(name = "Identity:  ",limits=c(60,100),
    colours=c("red","yellow","green","lightblue","darkblue"),
    breaks=c(60,80,100),labels=c("low","medium","high"),
    guide = guide_colorbar(title.theme=element_text(size=14, angle=0),title.vjust=1,
    barheight=0.6, barwidth=6, label.theme=element_text(size=10, angle=0)))
}
Cairo(width = 1600, height = 1600, 
      file="test7.pdf", type="pdf", pointsize=24, 
      bg = "transparent", canvas = "white", units = "px", dpi = 55)
do.call(grid.arrange,  plots_tqcover_cdna)
dev.off()

plots_tqcover_genome = list()
for(i in 1:length(dd$accession)){
  acc = dd$accession[i]
  dataset = assemblys_filter_genome[[acc]]
  setnames(dataset, gsub("%","",names(dataset)))
  dataset = filter(dataset,qLength>500)
  plots_tqcover_genome[[i]] = ggplot(data = dataset, aes(x = tCoverage, y = qCoverage, colour=identity)) + 
    geom_jitter(alpha=0.5) + geom_density2d() + theme(legend.position="bottom") +
    ggtitle(paste("tqCoverage (q>500)",paste(dd$accession[i],dd$tissue[i]))) +
    scale_colour_gradientn(name = "Identity:  ",limits=c(60,100),
    colours=c("red","yellow","green","lightblue","darkblue"),
    breaks=c(60,80,100),labels=c("low","medium","high"),
    guide = guide_colorbar(title.theme=element_text(size=14, angle=0),title.vjust=1,
    barheight=0.6, barwidth=6, label.theme=element_text(size=10, angle=0)))
}
Cairo(width = 1600, height = 1600, 
      file="before_75cutoff_tcoverage.pdf", type="pdf", pointsize=24, 
      bg = "transparent", canvas = "white", units = "px", dpi = 55)
do.call(grid.arrange,  plots_tqcover_genome)
dev.off()

