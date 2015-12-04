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
#mkdir -p $SEARCHPATH/filtersab
# rm $SEARCHPATH/filters/*
#
# 0.3. make data and llists
#
#echo "  making raw data and llist..." >&2
#cat ${SEARCHPATH}/transcripts_fpkm_3.fa | sed "s/>/>${SEARCHPATH}_/g" >${SEARCHPATH}/filters/raw_contigs.fa
#makellist ${SEARCHPATH}/filters/raw_contigs.fa >${SEARCHPATH}/filters/raw_contigs.llist
#
# 0.4. run blat using the CDS data
#
#echo "  running BLAT on ENSEMBL CDS..." >&2
#blat ${reference_cds} ${SEARCHPATH}/filters/raw_contigs.fa ${SEARCHPATH}/filters/blat_raw_cds.psl
#cat ${SEARCHPATH}/filters/blat_raw_cds.psl | ${refinery_blatparser} > ${SEARCHPATH}/filters/blat_raw_cds.best.tsv
#
# 0.5. re-filter by the CDS template
#
#echo "  filtering raw contigs using BLAT CDS alignment..." >&2
#Rscript ${refinery_bin}/filter_75_square.R ${SEARCHPATH}/filters/raw_contigs.llist ${SEARCHPATH}/filters/blat_raw_cds.best.tsv ${SEARCHPATH}/filters/filter0_cds_keepers.list
#include_mf ${SEARCHPATH}/filters/raw_contigs.fa ${SEARCHPATH}/filters/contigs_after_cds.fa ${SEARCHPATH}/filters/filter0_cds_keepers.list
#makellist ${SEARCHPATH}/filters/contigs_after_cds.fa >${SEARCHPATH}/filters/contigs_after_cds.llist
#
# 0.4. run blat using the CDNA data
#
#echo "  running BLAT on ENSEMBL CDNA..." >&2
#blat ${reference_cdna} ${SEARCHPATH}/filters/contigs_after_cds.fa ${SEARCHPATH}/filters/blat_cds_cdna.psl
#cat ${SEARCHPATH}/filters/blat_cds_cdna.psl | ${refinery_blatparser} > ${SEARCHPATH}/filters/blat_cds_cdna.best.tsv
#
# 0.5. re-filter by the CDNA template
#
#echo "  filtering raw contigs using BLAT CDNA alignment..." >&2
#Rscript ${refinery_bin}/filter_75_square.R ${SEARCHPATH}/filters/contigs_after_cds.llist ${SEARCHPATH}/filters/blat_cds_cdna.best.tsv ${SEARCHPATH}/filters/filter1_cdna_keepers.list
#include_mf ${SEARCHPATH}/filters/raw_contigs.fa ${SEARCHPATH}/filters/contigs_after_cdna.fa ${SEARCHPATH}/filters/filter1_cdna_keepers.list
#makellist ${SEARCHPATH}/filters/contigs_after_cdna.fa >${SEARCHPATH}/filters/contigs_after_cdna.llist
#
# 0.6 filter those with less than 20% overlap on REFSEQ PEP
#
#echo "  running BLASTX REFSEQ PEP alignment..." >&2
#blastall -a 3 -e .01 -d ${reference_pep} -p blastx -i ${SEARCHPATH}/filters/contigs_after_cdna.fa | gzip -9 - > ${SEARCHPATH}/filters/contigs_after_refseq_pep.blastx.out.gz
#zcat ${SEARCHPATH}/filters/contigs_after_refseq_pep.blastx.out.gz | ${refinery_bin}/blast_to_table.pl | ${refinery_bin}/hit_table_sorter.pl > ${SEARCHPATH}/filters/contigs_after_refseq_pep.best.tsv
#blat -t=dnax ${SEARCHPATH}/filters/contigs_after_cdna.fa -q=prot ${reference_pep_sequence} ${SEARCHPATH}/filters/blat_cdna_refseq_pep.psl
#cat ${SEARCHPATH}/filters/blat_cdna_refseq_pep.psl | ${refinery_blatparser} -f t > ${SEARCHPATH}/filters/blat_cdna_refseq_pep.best.tsv
#
#echo "  filtering contigs using BLAT REFSEQ PEP alignment..." >&2
#Rscript ${refinery_bin}/filter20_blast.R ${SEARCHPATH}/filters/contigs_after_cdna.llist ${SEARCHPATH}/filters/contigs_after_cdna_PEP.best.tsv ${SEARCHPATH}/filters/pep_keepers.list
#include_mf ${SEARCHPATH}/filters/raw_contigs.fa ${SEARCHPATH}/filters/contigs_after_pep.fa ${SEARCHPATH}/filters/pep_keepers.list
#makellist ${SEARCHPATH}/filters/contigs_after_pep.fa >${SEARCHPATH}/filters/contigs_after_pep.llist
#Rscript ${refinery_bin}/filter_75_square_inverse.R ${SEARCHPATH}/filters/contigs_after_cdna.llist ${SEARCHPATH}/filters/blat_cdna_refseq_pep.best.tsv ${SEARCHPATH}/filters/filter2_refseq_pep_keepers.list
#include_mf ${SEARCHPATH}/filters/raw_contigs.fa ${SEARCHPATH}/filters/contigs_after_refseq_pep.fa ${SEARCHPATH}/filters/filter2_refseq_pep_keepers.list
#makellist ${SEARCHPATH}/filters/contigs_after_refseq_pep.fa >${SEARCHPATH}/filters/contigs_after_refseq_pep.llist
#
#echo "  filtering contigs using BLAT REFSEQ PEP alignment..." >&2
#blat -t=dna ${reference_dna_sequence} -q=dna ${SEARCHPATH}/filters/contigs_after_refseq_pep.fa ${SEARCHPATH}/filters/blat_pep_refseq_dna.psl
#cat ${SEARCHPATH}/filters/blat_pep_refseq_dna.psl | ${refinery_blatparser} > ${SEARCHPATH}/filters/blat_pep_refseq_dna.best.tsv
#Rscript ${refinery_bin}/filter_75_square.R ${SEARCHPATH}/filters/contigs_after_refseq_pep.fa ${SEARCHPATH}/filters/blat_pep_refseq_dna.best.tsv ${SEARCHPATH}/filters/filter3_refseq_dna_keepers.list
#include_mf ${SEARCHPATH}/filters/raw_contigs.fa ${SEARCHPATH}/filters/contigs_after_refseq_dna.fa ${SEARCHPATH}/filters/filter3_refseq_dna_keepers.list
#makellist ${SEARCHPATH}/filters/contigs_after_refseq_dna.fa >${SEARCHPATH}/filters/contigs_after_refseq_dna.llist
#
echo "  filtering contigs using BLAT GENOME alignment..." >&2
blat -t=dna ${reference_genome} -q=dna ${SEARCHPATH}/filters/contigs_after_refseq_pep.fa ${SEARCHPATH}/filters/blat_pep_genome.psl
cat ${SEARCHPATH}/filters/blat_pep_genome.psl | ${refinery_blatparser} > ${SEARCHPATH}/filters/blat_pep_genome.best.tsv
Rscript ${refinery_bin}/filter_75_square.R ${SEARCHPATH}/filters/contigs_after_refseq_pep.fa ${SEARCHPATH}/filters/blat_pep_genome.best.tsv ${SEARCHPATH}/filters/filter4_genome_keepers.list
include_mf ${SEARCHPATH}/filters/raw_contigs.fa ${SEARCHPATH}/filters/contigs_after_genome.fa ${SEARCHPATH}/filters/filter4_genome_keepers.list
makellist ${SEARCHPATH}/filters/contigs_after_genome.fa >${SEARCHPATH}/filters/contigs_after_genome.llist
