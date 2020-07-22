#!/bin/bash

if [ "$#" -lt 4 ]; then
    echo "Usage: $0 buildDir runName manVersion[B4,C0] classVar[Sample_Class,Sample_Name]"
    exit 1
fi
buildDir=$1
runName=$2
version=$3
classVar=$4

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
prgmDir="covic"
prgmTag="merge_builds"

TOP_MAC=/Users/bbarnes/Documents/CustomerFacing
TOP_LIX=/illumina/scratch/darkmatter/Projects/COVIC

if [ -e ${TOP_MAC} ]; then
    TOP=${TOP_MAC}
    SRC=${TOP_MAC}/${prgmTop}
    DAT=${SRC}/dat
    CONDA=mac
    RSCRIPT=/usr/local/bin/Rscript

    sampleCsv=${TOP}/sampleSheets/annotation/Human-Classification_COVID_Count-656_AnnotatedMultiSampleSheet.csv

elif [ -e ${TOP_LIX} ]; then
    TOP=${TOP_LIX}
    SRC=${TOP_LIX}/tools/${prgmTop}
    DAT=${SRC}/dat

    CONDA=conda_4.6.8
    # CONDA=Anaconda2-2019.10-Linux-x86_64
    # CONDA=Anaconda3-2019.10-Linux-x86_64
    RSCRIPT=/illumina/scratch/darkmatter/thirdparty/${CONDA}/bin/Rscript

    sampleCsv=${TOP}/sampleSheets/annotation/Human-Classification_COVID_Count-656_AnnotatedMultiSampleSheet.csv

else
    echo "Unrecognized top directory!"
    exit
fi

EXE=${SRC}/scripts/R/${prgmDir}/${prgmTag}.R

CMD=${RSCRIPT}" "${EXE}
CMD+=" --"Rscript=${RSCRIPT}

# Directories::
outDir=${TOP}/${prgmTag}

CMD+=" --buildDir"=${buildDir}
CMD+=" --runName"=${runName}
CMD+=" --outDir"=${outDir}

platform='EPIC'
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
