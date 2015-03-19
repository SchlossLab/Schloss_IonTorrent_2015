#!/bin/bash

################################################################################
#
# basic_error_analysis.sh
#
#
# Here we get the error rates for the individual reads with no curation. This is
# similar to what we did with the MiSeq analysis:
#   https://github.com/SchlossLab/Kozich_MiSeqSOP_AEM_2013/blob/master/code/single_read_analysis.sh
#
# Dependencies...
# * The *.fasta and *.qual files out of the data/enzyme?/ directory
#
# Produces
# * The error rates for each individual read including
#   - *.error.summary
#   - *.error.matrix
#   - *.error.quality
#   - *.error.rev.qual
#   - *.error.rev.seq
#   - *.trim.summary
# * The starting and ending coordinates for every aligned read
# * These will be used to make box plots showing the relationship between the
#   quality scores and error type, the effect of the position in the read on
#   the error rate, and the substitution matrix
#
################################################################################

FASTA=$1
FASTA=data/raw_enzyme1/soil1.fasta
QUAL=$(echo $FASTA | sed -e s/fasta/qual/)

RAW_PATH=$(echo $FASTA | sed -e 's/\/\w*.fasta//')
BASIC_PATH=$(echo $RAW_PATH | sed -e 's/raw/basic/')

FILE_STUB=$(echo $FASTA | sed -e 's/.*\/\(\w*\).fasta/\1/')
BASIC_STUB=$BASIC_PATH/$FILE_STUB

REF=$BASIC_STUB.HMP_MOCK.align
FILTER_REF=$(echo $REF | sed -e s/align/filter.fasta/)

REPORT=$BASIC_STUB.trim.align.report

cp data/references/HMP_MOCK.v35.align $REF

mothur "#set.dir(output=$BASIC_PATH);
    trim.seqs(fasta=$FASTA, qfile=$QUAL, oligos=data/iontorrent.oligos, bdiffs=1, pdiffs=2, flip=T, processors=8);
    align.seqs(fasta=current, reference=$REF);
    summary.seqs();
    filter.seqs(fasta=current-$REF, vertical=T);
    seq.error(fasta=current, qfile=current, report=$REPORT, reference=$FILTER_REF);"

rm $BASIC_STUB.trim.fasta
rm $BASIC_STUB.trim.qual
rm $BASIC_STUB.scrap.fasta
rm $BASIC_STUB.scrap.qual
rm $BASIC_STUB.groups
rm $BASIC_STUB.trim.align
rm $BASIC_STUB.trim.align.report
rm $BASIC_STUB.trim.flip.accnos
rm $BASIC_STUB*.filter
rm $BASIC_STUB.trim.filter.fasta
rm $BASIC_STUB.HMP_MOCK.filter.fasta
rm $BASIC_STUB.trim.filter.error.seq
rm $BASIC_STUB.trim.filter.error.chimera
rm $BASIC_STUB.trim.filter.error.qual.forward
rm $BASIC_STUB.trim.filter.error.seq.forward
rm $BASIC_STUB.trim.filter.error.count
rm $BASIC_STUB.trim.filter.error.ref
rm $BASIC_STUB.HMP_MOCK.8mer
rm $BASIC_STUB.HMP_MOCK.align
