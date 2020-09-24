#!/bin/bash

if [ "$#" -lt 5 ]; then
    echo "Usage: $0 outTopDir buildDir runName manVersion[B4,C0] classVar[Sample_Class,Sample_Name]"
    exit 1
fi
outTop=$1
buildDir=$2
runName=$3
version=$4
classVar=$5

# version='B4'
# version='C0'
# classVar="Sample_Class"

# Parallel/Cluster Options::
# single=true
single=false
parallel=true
cluster=true

verbose=3

#
# Standard Options:: DON'T CHANGE::
#
prgmTop="Infinium_Methylation_Workhorse"
prgmDir="analysis"
prgmTag="merge_builds"

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

SRC=${TOP}/${prgmTop}
DAT=${SRC}/dat
EXE=${SRC}/scripts/R/${prgmDir}/${prgmTag}.R

# sampleCsv=${TOP}/sampleSheets/annotation/Human-Classification_COVID_Count-656_AnnotatedMultiSampleSheet.csv
# sampleCsv=${DAT}/sampleSheets/COVIC/Human-Classification_COVID_Count-921_AnnotatedOnlySampleSheet.csv.gz

CMD=${RSCRIPT}" "${EXE}
CMD+=" --"Rscript=${RSCRIPT}

# Directories::
outDir=${outTop}/${prgmTag}

CMD+=" --buildDir"=${buildDir}
CMD+=" --runName"=${runName}
CMD+=" --outDir"=${outDir}

platform='EPIC'
CMD+=" --platform"=${platform}
CMD+=" --version"=${version}

CMD+=" --classVar"=${classVar}
# CMD+=" --sampleCsv"=${sampleCsv}

if [ "${single}" = true ]; then
    CMD+=" --single"
fi
if [ "${parallel}" = true ]; then
    CMD+=" --parallel"
fi
if [ "${cluster}" = true ]; then
    CMD+=" --cluster"
fi

# CMD+=" --"=${}

# Verbosity Options::
CMD+=" --"verbose=${verbose}

# mkdir -p ${outDir}
echo ${CMD}
echo

${CMD}

## End of file
