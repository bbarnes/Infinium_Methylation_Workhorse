#!/bin/bash

if [ "$#" -lt 7 ]; then
    echo "Usage: $0 mergeDir(s) modelDir runName classVar[Sample_Class,Sample_Name] trainClass[nSARSCov2,pSARSCov2] version[B4,C0] seedStr[seed-13,seed-42,seed-61]"
    exit 1
fi
mergeDir=$1
modelDir=$2
runName=$3
classVar=$4
trainClass=$5
version=$6
seedStr=$7

# General Parameters::
#
# classVar="Sample_Class"
# trainClass="nSARSCov2,pSARSCov2"
# classVar="Sample_Name"
# trainClass="HELA,JURKAT,MCF7,nSARSCov2,RAJI"
# classVar="Karyotype_1_Call"
# trainClass="XaXi,XaY"

platform="EPIC"

lociBetaKey="i_beta,ind_beta"
lociPvalKey="i_poob"
lociPvalMin="1.0,0.1"

# Parallel/Cluster Options::
# single=true
single=false
parallel=true
cluster=true

verbose=3

#
# Standard Options:: DON'T CHANGE:: 
#
prgmTop="workhorse"
prgmDir="covic"
prgmTag="predict_builds"

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
outDir=${TOP}/${prgmTag}/${platform}/${version}

# idatsDir=NULL
CMD+=" --mergeDir"=${mergeDir}
CMD+=" --modelDir"=${modelDir}
CMD+=" --runName"=${runName}
CMD+=" --outDir"=${outDir}

CMD+=" --seed_dir"=${seedStr}

CMD+=" --classVar"=${classVar}
CMD+=" --trainClass"=${trainClass}

CMD+=" --lociBetaKey"=${lociBetaKey}
CMD+=" --lociPvalKey"=${lociPvalKey}
CMD+=" --lociPvalKey"=${lociPvalKey}
CMD+=" --lociPvalMin"=${lociPvalMin}

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
