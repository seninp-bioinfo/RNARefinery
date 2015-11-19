#!/bin/sh
# Sample shell script to perform filtering
echo "Reading refinery config..." >&2
source refinery.properties
#
# parse command line args
#
while [[ $# > 1 ]]
do
key="$1"
case $key in
    -s|--searchpath)
    SEARCHPATH="$2"
    shift # past argument
    ;;
    --default)
    DEFAULT=YES
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done
echo SEARCH PATH USED "$SEARCHPATH" >&2
#
# 0.1. making sure FPKM3 exists
#
if [ ! -f $SEARCHPATH/transcripts_fpkm_3.fa ]; then
    echo "File not found, breaking!"
    exit 10
fi
echo "  $SEARCHPATH/transcripts_fpkm_3.fa found..." >&2
#
# 0.2. make the filter folder
#
mkdir -p $SEARCHPATH/filters
#
# 0.3. make data and llists
#
echo "  making raw data and llist..." >&2
cat ${SEARCHPATH}/transcripts_fpkm_3.fa | sed "s/>/>${SEARCHPATH}_/g" >${SEARCHPATH}/filters/raw_contigs.fa
makellist ${SEARCHPATH}/filters/raw_contigs.fa >${SEARCHPATH}/filters/raw_contigs.llist
#
# 0.3. run blat using the data
#
echo "  running BLAT on ENSEMBL CDS..." >&2
blat ${reference_cds} ${SEARCHPATH}/filters/raw_contigs.fa ${SEARCHPATH}/filters/blat_raw_cds.psl
cat ${SEARCHPATH}/filters/blat_raw_cds.psl | blat_parser.pl > ${SEARCHPATH}/filters/blat_raw_cds.best.tsv
#
# 0.4 filter those with less than 20% overlap on CDS
#
echo "  filtering raw contigs using BLAT CDS alignment..." >&2
Rscript ${refinery_bin}/filter20_blat.R ${SEARCHPATH}/filters/raw_contigs.llist ${SEARCHPATH}/filters/blat_raw_cds.best.tsv ${SEARCHPATH}/filters/cds_keepers.list
include_mf ${SEARCHPATH}/filters/raw_contigs.fa ${SEARCHPATH}/filters/contigs_after_cds.fa ${SEARCHPATH}/filters/cds_keepers.list
makellist ${SEARCHPATH}/filters/contigs_after_cds.fa >${SEARCHPATH}/filters/contigs_after_cds.llist
#
# 0.5 filter those with less than 20% overlap on CDNA
#
echo "  running BLAT on ENSEMBL CDNA..." >&2
blat ${reference_cdna} ${SEARCHPATH}/filters/contigs_after_cds.fa ${SEARCHPATH}/filters/blat_cds_cdna.psl
cat ${SEARCHPATH}/filters/blat_cds_cdna.psl | blat_parser.pl > ${SEARCHPATH}/filters/blat_cds_cdna.best.tsv
#
echo "  filtering raw contigs using BLAT CDNA alignment..." >&2
Rscript ${refinery_bin}/filter20_blat.R ${SEARCHPATH}/filters/contigs_after_cds.llist ${SEARCHPATH}/filters/blat_cds_cdna.best.tsv ${SEARCHPATH}/filters/cdna_keepers.list
include_mf ${SEARCHPATH}/filters/raw_contigs.fa ${SEARCHPATH}/filters/contigs_after_cdna.fa ${SEARCHPATH}/filters/cdna_keepers.list
makellist ${SEARCHPATH}/filters/contigs_after_cdna.fa >${SEARCHPATH}/filters/contigs_after_cdna.llist
#
# 0.5 filter those with less than 20% overlap on REFSEQ PEP
#
blastall -a 5 -e .01 -d ${reference_pep} -p blastx -i ${SEARCHPATH}/filters/contigs_after_cdna.fa | gzip -9 - > ${SEARCHPATH}/filters/contigs_after_cdna_PEP.blastx.out.gz
zcat ${SEARCHPATH}/filters/contigs_after_cdna_PEP.blastx.out.gz | ${refinery_bin}/blast_to_table.pl | ${refinery_bin}/hit_table_sorter.pl > ${SEARCHPATH}/filters/contigs_after_cdna_PEP.best.tsv

