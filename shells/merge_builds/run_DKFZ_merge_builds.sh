#!/bin/bash

if [ "$#" -lt 6 ]; then
    echo "Usage: $0 outTopDir buildDir runName platform[HM450,EPIC] manVersion[B4,C0] classVar[Sample_Class,Sample_Name]"
    exit 1
fi
outTop=$1
buildDir=$2
runName=$3
platform=$4
version=$5
classVar=$6

# platform='EPIC'
# version='B4'
# version='C0'
# classVar="Sample_Class"

# Parallel/Cluster Options::
# single=true
single=false
parallel=true
cluster=true

verbose=30

#
# Standard Options:: DON'T CHANGE::
#
prgmTop="Infinium_Methylation_Workhorse"
prgmDir="analysis"
prgmTag="merge_builds"

TOP_MAC1=/Users/bbarnes/Documents/Projects/methylation/tools
TOP_MAC2=/Users/bretbarnes/Documents/tools
TOP_LIX1=/illumina/scratch/darkmatter/tools
TOP_LIX2=/illumina/scratch/darkmatter/tools/latest

RSCRIPT=/usr/local/bin/Rscript

if [ -e ${TOP_MAC1} ]; then
    TOP=${TOP_MAC1}

elif [ -e ${TOP_MAC2} ]; then
    TOP=${TOP_MAC2}

elif [ -e ${TOP_LIX1} ]; then
    TOP=${TOP_LIX1}

    CONDA=conda_4.6.8
    # CONDA=Anaconda2-2019.10-Linux-x86_64
    # CONDA=Anaconda3-2019.10-Linux-x86_64
    RSCRIPT=/illumina/scratch/darkmatter/thirdparty/${CONDA}/bin/Rscript

elif [ -e ${TOP_LIX2} ]; then
    TOP=${TOP_LIX2}

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
sampleCsv=${DAT}/sampleSheets/DKFZ/dkfz_sentrix-name_classes.csv.gz

CMD=${RSCRIPT}" "${EXE}
CMD+=" --"Rscript=${RSCRIPT}

# Directories::
outDir=${outTop}/${prgmTag}

CMD+=" --buildDir"=${buildDir}
CMD+=" --runName"=${runName}
CMD+=" --outDir"=${outDir}

CMD+=" --platform"=${platform}
CMD+=" --version"=${version}

CMD+=" --classVar"=${classVar}
CMD+=" --sampleCsv"=${sampleCsv}

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