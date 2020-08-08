#!/bin/bash

if [ "$#" -lt 9 ]; then
    echo "Usage: $0 fasta outDir runName platform version build genome des1_csv des2_csv MAX"
    exit 1
fi
fasta=$1
outDir=$2
runName=$3

platform=$4
version=$5
build=$6
genome=$7

des1_csv=$8
des2_csv=$9
max=${10}

verbose=3

# Program Variables::
#
prgmTop="Infinium_Methylation_Workhorse"
prgmDir="probe_design"
prgmTag="tile_main"

TOP_MAC=/Users/bbarnes/Documents/Projects/methylation/tools
TOP_LIX=/illumina/scratch/darkmatter/Projects/COVIC/tools

if [ -e ${TOP_MAC} ]; then
    TOP=${TOP_MAC}
    CONDA=mac
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

SRC=${TOP}/${prgmTop}
EXE=${SRC}/scripts/R/${prgmDir}/${prgmTag}.R

CMD=${RSCRIPT}" "${EXE}
CMD+=" --"Rscript=${RSCRIPT}

CMD+=" --outDir"=${outDir}
CMD+=" --runName"=${runName}
CMD+=" --fasta"=${fasta}
CMD+=" --platform"=${platform}
CMD+=" --version"=${version}
CMD+=" --build"=${build}
CMD+=" --genome"=${genome}

CMD+=" --des1_csv"=${des1_csv}
CMD+=" --des2_csv"=${des2_csv}

CMD+=" --max"=${max}

# Verbosity Options::
CMD+=" --"verbose=${verbose}

# mkdir -p ${outDir}
echo ${CMD}

${CMD}

## End of file