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
# 0.3. run blat using the data, re-filter by the template
#
echo "  filtering raw contigs using BLAT CDS alignment..." >&2
Rscript ${refinery_bin}/filter80_tCover.R ${SEARCHPATH}/filters/raw_contigs.llist ${SEARCHPATH}/filters/blat_raw_cds.best.tsv ${SEARCHPATH}/filters/tcover_cds_keepers.list
include_mf ${SEARCHPATH}/filters/raw_contigs.fa ${SEARCHPATH}/filters/tcover_contigs_after_cds.fa ${SEARCHPATH}/filters/tcover_cds_keepers.list
makellist ${SEARCHPATH}/filters/tcover_contigs_after_cds.fa >${SEARCHPATH}/filters/tcover_contigs_after_cds.llist