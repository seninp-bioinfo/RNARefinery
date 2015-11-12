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
require(seqinr)
#
properties = read.table("/home/psenin/.rnarefinery/db.properties", header = F, as.is = T)
db_credentials = data.frame(t(properties$V2), stringsAsFactors = F); names(db_credentials) = properties$V1
session <- dbConnect(MySQL(), host = db_credentials$host, dbname = db_credentials$db, 
                     user = db_credentials$user, password = db_credentials$password)
#
#tissue = "embryo"
#tissue = "breast muscle"
tissue = "tibial bone"
runs = dbGetQuery(session, paste("select * from sra_runs where UPPER(tissue)=",shQuote(toupper(tissue)),
                  " and contigs>0 and processed=0"))
#
#out_fname = "raw_contigs_embryo.fa"
#out_fname = "raw_contigs_bmuscle.fa"
out_fname = "raw_contigs_tbone.fa"
#
for (id in runs$id) {
  contigs = dbGetQuery(session, paste("select * from raw_contigs where run_id=", id))
  clist = dlply(contigs, .(name), function(x) {
    as.SeqFastadna(s2c(x$sequence), name = x$name, Annot = NULL)
  })
  write.fasta(clist, names = attributes(clist)$name, nbchar = 80, file.out = out_fname, open = "a")
}
