#
# if runs accession in the list this updates the date
#
require(RMySQL)
#
require(reshape)
require(plyr)
require(dplyr)
require(stringr)
#
properties = read.table("/home/psenin/.rnarefinery/db.properties", header = F, as.is = T)
db_credentials = data.frame(t(properties$V2), stringsAsFactors = F); names(db_credentials) = properties$V1
session <- dbConnect(MySQL(), host = db_credentials$host, dbname = db_credentials$db, 
                     user = db_credentials$user, password = db_credentials$password)
#
runs = read.table("/work2/project/sigenae/Project_RNA_Refinery/pipelineOut/00_processed_runs.list", header = F, as.is = T)
#
a_ply(runs, 1, function(x) {
  dbGetQuery(session, paste("update sra_runs set processed=",shQuote("2015-09-01 00:00:00"),
                            " where run_accession=", shQuote(x)))
})
