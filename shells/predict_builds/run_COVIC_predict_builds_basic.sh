#!/bin/bash

if [ "$#" -lt 4 ]; then
    echo "Usage: $0 dataDir runName targetName outDir"
    exit 1
fi
DAT=$1
RUN=$2
TAR=$3
OUT=$4

CHIP="EPIC"
VER="C0"
CLASS="Sample_Class"
trainClass="nSARSCov2,pSARSCov2"

BETA_KEY="ind_beta"
PVAL_KEY="i_poob"
PVAL_MIN="1.0,0.1"

VERBOSE=30

# RUN=COVIC-Set-15
# RUN=COVIC-Set-17
# RUN=COVIC-Set-57
# RUN=COVIC-Set5-10062020
# RUN=COVIC-Set1-15052020
# RUN=COVIC-Set7-06082020

# RUN=COVIC-Set7-06082020
# TAR=COVIC-Set-15

prgmTop="Infinium_Methylation_Workhorse"
prgmDir="analysis"
prgmTag="predict_models"

TOP_MAC=/Users/bbarnes/Documents/Projects/methylation/tools
TOP_LIX=/illumina/scratch/darkmatter/tools

if [ -e ${TOP_MAC} ]; then
    TOP=${TOP_MAC}
    RSCRIPT=/usr/local/bin/Rscript

elif [ -e ${TOP_LIX} ]; then
    TOP=${TOP_LIX}

    CONDA=conda_4.6.8
    # CONDA=Anaconda2-2019.10-Linux-x86_64
    # CONDA=Anaconda3-2019.10-Linux-x86_64
    RSCRIPT=/illumina/scratch/darkmatter/thirdparty/${CONDA}/bin/Rscript

else
    echo "Unrecognized top directory!"
    exit
fi

EXE=${SRC}/scripts/R/${prgmDir}/${prgmTag}.R
MOD=${DAT}/build_models/${CHIP}/${VER}/${CLASS}/${RUN}
MER=${DAT}/merge_builds/${CHIP}/${VER}/${CLASS}/${TAR}

echo "EXE="${EXE}
echo "DAT="${DAT}
echo "MOD="${MOD}
echo "MER="${MER}
echo "OUT="${OUT}

${RSCRIPT} ${EXE} -o ${OUT} -m ${MER} -r ${RUN} \
    --modelDir=${MOD} \
    --classVar=${CLASS} \
    --trainClass=${TRAIN} \
    --lociBetaKey=${BETA_KEY} \
    --lociPvalKey=${PVAL_KEY} \
    --lociPvalMin=${PVAL_MIN} \
    --Rscript=${RSCRIPT} \
    --verbose=${VERBOSE}


# End of file
