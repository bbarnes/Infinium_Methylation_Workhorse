

TOP="/Users/bretbarnes/Documents/data/improbe/scratch/cgnDB/dbSNP_Core4/design-input"

BLD="GRCm10"
BLD="GRCh37"

DIR="${TOP}/min"
INP="${TOP}/${BLD}.cgn.bed.gz"
OUT="${DIR}/${BLD}.chr-pos-srd.slim.cgn-sorted.txt.gz"

mkdir -p ${DIR}

gzip -dc ${INP} | \
    cut -f 1,3,4,6 | \
    perl -pe 's/^chr//' | \
    sort -n -k 3,3 | \
    gzip -c -> ${OUT}


echo "done."

exit

#    head | \

# End of file
