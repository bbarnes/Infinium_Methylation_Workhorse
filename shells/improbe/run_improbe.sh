#!/bin/bash

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 idatsDir"
    exit 1
fi

INP=$1
DES=${INP}.improbe-design.tsv
LOG=${INP}.improbe-design.log

TOP=/illumina/scratch/darkmatter
DAT=${TOP}/dat
EXE=${TOP}/bin/improbe

TAN=${DAT}/Tango_A_or_B_11mer_s1.dat
MER=${DAT}/human-36.1-methyl-and-unmethyl-13mer-s3-for-infinium-methylation.dat

echo ${EXE} -oASPE -tBOTH -cBoth -n${MER} -a${TAN} -V ${INP} > ${DES} 2> ${LOG}

${EXE} -oASPE -tBOTH -cBoth -n${MER} -a${TAN} -V ${INP} > ${DES} 2> ${LOG}

echo "Done."

# End of file
