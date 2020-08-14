#!/bin/bash

if [ "$#" -lt 8 ]; then
    echo "Usage: $0 mergeDir(s) sampleSheet runName version[C0,B4] classVar[Sample_Class,Sample_Name] trainClass[nSARSCov2,pSARSCov2] lociBetaKey[i_beta,ind_beta] localPvalKey[i_poob,i_negs] lociPvalMin[1.0, 0.2, 0.1, 0.01] featureDMLs[100,1000] seeds[13,42]"
    exit 1
fi

mergeDir=$1

runName=$2
version=$3
classVar=$4
trainClass=$5

lociBetaKey=$6
lociPvalKey=$7
lociPvalMin=$8
featuresDml=$9
seeds=${10}

platform="EPIC"

# General Parameters::
#
# classVar="Sample_Class"
# trainClass="nSARSCov2,pSARSCov2"

samplePvalName="Poob_Pass_0_Perc"
samplePvalPerc=96

# lociBetaKey="i_beta,ind_beta"
# lociPvalKey="i_poob"
# lociPvalMin="1.0,0.1"
# lociPvalMin="0.1"

# featuresDml="100,1000,10000,15000"
# featuresDml="20000"
# featuresDml="100,1000"
featuresDbl=100

# seeds="13,17,41,42,43,61,69"
# seeds="13,42"

# Parallel/Cluster Options::
# single=true
single=false
# parallel=true
parallel=false
cluster=true

verbose=3

#
# Standard Options:: DON'T CHANGE::
#
prgmTop="Infinium_Methylation_Workhorse"
prgmDir="analysis"
prgmTag="build_models"

TOP_MAC=/Users/bbarnes/Documents/Projects/methylation/tools
TOP_LIX=/illumina/scratch/darkmatter/Projects/COVIC/tools

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

featuresCsv=${DAT}/sampleSheets/dmls/Ivana-145.csv.gz

CMD=${RSCRIPT}" "${EXE}
CMD+=" --"Rscript=${RSCRIPT}

# Directories::
outDir=${TOP}/${prgmTag}/${platform}/${version}

# idatsDir=NULL
CMD+=" --mergeDir"=${mergeDir}
CMD+=" --runName"=${runName}
CMD+=" --outDir"=${outDir}

# platform='EPIC'
# CMD+=" --platform"=${platform}
# CMD+=" --version"=${version}

CMD+=" --classVar"=${classVar}
CMD+=" --trainClass"=${trainClass}

CMD+=" --samplePvalName"=${samplePvalName}
CMD+=" --samplePvalPerc"=${samplePvalPerc}

CMD+=" --lociBetaKey"=${lociBetaKey}
CMD+=" --lociPvalKey"=${lociPvalKey}
CMD+=" --lociPvalKey"=${lociPvalKey}
CMD+=" --lociPvalMin"=${lociPvalMin}

CMD+=" --featuresCsv"=${featuresCsv}
CMD+=" --featuresDml"=${featuresDml}
CMD+=" --featuresDbl"=${featuresDbl}
CMD+=" --seeds"=${seeds}

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