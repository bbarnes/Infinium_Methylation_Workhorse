#!/bin/bash

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 genomeBuild"
    exit 1
fi
BUILD=$1

# TOP_DIR="/Users/bretbarnes/Documents"
# DAT_DIR="${TOP_DIR}/data/improbe/GRCh36-38.mm10.FWD_SEQ.21022021/designOutput"

TOP_DIR="/illumina/scratch/darkmatter/Projects/dbCGN"
DAT_DIR="${TOP_DIR}/data/cgnDB/GRCh36-38.mm10.FWD_SEQ.21022021/designOutput"

OUT_DIR="${DAT_DIR}/prb49"

DES="${DAT_DIR}/${BUILD}-21022021_improbe-designOutput.tsv.gz"
OUT_U50="${OUT_DIR}/${BUILD}-21022021_improbe-designOutput.m50-cgn.seq-sorted.tsv.gz"
OUT_U49="${OUT_DIR}/${BUILD}-21022021_improbe-designOutput.m49-cgn.seq-sorted.tsv.gz"

echo "BUILD=${BUILD}"
echo "DAT_DIR=${DAT_DIR}"
echo "OUT_DIR=${OUT_DIR}"
echo ""
echo "DES=${DES}"
echo "OUT_U50=${OUT_U50}"
echo "OUT_U49=${OUT_U49}"
echo ""

mkdir -p ${OUT_DIR}

gzip -dc ${DES} | cut -f 1,21-23,16 | \
    # head -n 20 | \
    grep -v Seq_ID | \
    perl -pe 's/\n$//; @d=split(/\t/,$_); s/^.*$//; $d[3]=substr($d[3],0,1); $seq=substr($d[1],0,49); $e=substr($d[1],49);  print "$seq\t$e\t$d[0]\t$d[2]\t$d[3]\t$d[4]\n"' | \
    sort -u -k 1,3 | \
    gzip -c -> ${OUT_U50}

exit

gzip -dc ${OUT_U50} | \
    perl -pe 's/^[A-Z]//;' | \
    sort -u -k 1,2 | \
    gzip -c -> ${OUT_U49}

echo "done."

# End of file
