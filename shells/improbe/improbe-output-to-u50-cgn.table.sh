#!/bin/bash

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 genomeBuild"
    exit 1
fi
BUILD=$1

TOP_DIR="/Users/bretbarnes/Documents"
DAT_DIR="${TOP_DIR}/data/improbe/GRCh36-38.mm10.FWD_SEQ.21022021/designOutput"
OUT_DIR="${DAT_DIR}/prbDB"

DES="${DAT_DIR}/${BUILD}-21022021_improbe-designOutput.tsv.gz"
OUT="${OUT_DIR}/${BUILD}-21022021_improbe-designOutput.u50-cgn.seq-sorted.tsv.gz"

echo "BUILD=${BUILD}"
echo "DAT_DIR=${DAT_DIR}"
echo "OUT_DIR=${OUT_DIR}"
echo ""
echo "DES=${DES}"
echo "OUT=${OUT}"
echo ""

mkdir -p ${OUT_DIR}

gzip -dc ${DES} | cut -f 1,21-23,48 | \
    # head -n 20 | \
    grep -v Seq_ID | \
    perl -pe 's/\n$//; @d=split(/\t/,$_); s/^.*$//; $d[2]=substr($d[2],0,1);  print "$d[4]\t$d[0]\t$d[1]\t$d[2]\t$d[3]\n"' | \
    sort -u -k 1,2 | \
    gzip -c -> ${OUT}

echo "done."

# End of file
