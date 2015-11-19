## get the list of completed assemblys
#
setwd("/work2/project/sigenae/Project_RNA_Refinery/SRAData/")
datafiles <-
  data.frame(files = list.files(".", "*.fastq*", include.dirs = F))

library(stringr)
datafiles$run = str_extract(datafiles$files,"[^_]*")

library(plyr)
ddply(datafiles, .(run), function(x) {
  print(paste(x))
})

library(dplyr)
pre_run <-
  datafiles %>% group_by(run) %>% summarise_each(funs(first(.), last(.)))

setwd("/work2/project/sigenae/Project_RNA_Refinery/runDrapOut/")
assembled <- data.frame(runs = list.files(".", include.dirs = T))

pre_run <- pre_run[!(pre_run$run %in% assembled$runs),]

# /home/sigenae/bin/runDrap -1 /work2/project/sigenae/Project_RNA_Refinery/SRAData/
# SRR1811333_1.fastq.gz
# -2 /work2/project/sigenae/Project_RNA_Refinery/SRAData/
# SRR1811333_2.fastq.gz
# -o /work2/project/sigenae/Project_RNA_Refinery/runDrapOut/
# SRR1811333
# --alignTrim

pre_run$command = paste(
  "/home/sigenae/bin/runDrap -1 /work2/project/sigenae/Project_RNA_Refinery/SRAData/",
  pre_run$first,
  " -2 /work2/project/sigenae/Project_RNA_Refinery/SRAData/",
  pre_run$last,
  " -o /work2/project/sigenae/Project_RNA_Refinery/runDrapOut/",
  pre_run$run, " --alignTrim", sep = ""
)

setwd("/work2/project/sigenae/Project_RNA_Refinery/")

write.table(
  pre_run$command, "rundrap_2015_10_20.sh", col.names = F, row.names = F, quote = F
)