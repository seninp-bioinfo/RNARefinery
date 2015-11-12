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
db_credentials = as.data.frame(t(properties$V2)); names(db_credentials) = properties$V1
session <- dbConnect(MySQL(), host = db_credentials$hostdb_credentials, db = db_credentials$db, 
                     user = db_credentials$user, password = db_credentials$password)
#
runs = dbGetQuery(session, "select * from sra_runs")
#
a_ply(runs$run_accession, 1, function(x) {
  count = dbGetQuery(session, paste("select count(*) from raw_contigs where name like ", 
                                    shQuote(paste(x,"%",sep = ''))))
  if (count > 0) {
    dbGetQuery(session, paste("update sra_runs set contigs=",count," where run_accession=", shQuote(x)))
  }  
})
