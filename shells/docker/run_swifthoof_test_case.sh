#!/bin/bash

BRD=Infinium_Methylation_Workhorse_Centos
EXP=TestCase
VER=v.1.1
SRC="bbarnesimdocker/im_workhorse:${BRD}.${VER}"
CMD=run_swifthoof.sh

DAT=/repo/Infinium_Methylation_Workhorse/dat/idats_${EXP}
OUT=./docker-${BRD}-${VER}

VERBOSE=4
WORKFLOW="nd,ind"

# mkdir -p ${OUT}

# docker run -i --rm -v ${DAT}:/input -v ${OUT}:/output -w /work ${SRC} ${CMD} \
docker run -i --rm -v ${OUT}:/output -w /work ${SRC} ${CMD} \
       --runName=${EXP} \
       --workflow="${WORKFLOW}" \
       --write_call \
       --write_sigs \
       --load_sset \
       --load_idat \
       --save_sset \
       --save_idat \
       --pval="pOOBAH,PnegEcdf" --minPval="0.1,0.02" --minPerc="90,98" \
       --verbose=${VERBOSE}

echo "swifthoof done..."

# End of file
