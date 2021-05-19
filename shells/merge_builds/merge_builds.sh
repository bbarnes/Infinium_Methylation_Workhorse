#!/bin/bash

EXP="EPIC-8x1-EM-Sample-Prep"
EXP="NA12878.v0"

# TOP="/Users/bretbarnes/Documents"
# SRC="${TOP}/tools/Infinium_Methylation_Workhorse"
# EXE="${SRC}/scripts/R/

EXE_NAME="Infinium_Methylation_Workhorse/scripts/R/analysis/merge_builds.R"
EXE_A="/repo/${EXE_NAME}"
EXE_B="/Users/bretbarnes/Documents/tools/${EXE_NAME}"

RSCRIPT_A=/usr/local/bin/Rscript
RSCRIPT_B=/usr/bin/Rscript

if [ -e ${RSCRIPT_A} ]; then
    RSCRIPT=${RSCRIPT_A}
elif [ -e ${RSCRIPT_B} ]; then
    RSCRIPT=${RSCRIPT_B}
else
    echo "Unrecognized Rscript EXE!"
    exit
fi

echo "RSC="${RSCRIPT}

if [ -e ${EXE_A} ]; then
    EXE=${EXE_A}

    echo "Docker Run..."
    echo "EXE="${EXE}

    INP=/input
    OUT=/output
    CMD="${RSCRIPT} ${EXE} --Rscript ${RSCRIPT} --datDir=${INP} --outDir=${OUT} $@"
    
elif [ -e ${EXE_B} ]; then
    EXE=${EXE_B}
    
    echo "Local Run..."
    echo "EXE="${EXE}
    
    INP=$1
    OUT=$2

    if [ "$#" -lt 2 ]; then
	echo "Usage: $0 buildDir outDir [options]"
	exit 1
    fi
    
    CMD="${RSCRIPT} ${EXE} --Rscript ${RSCRIPT} $@"
    
else
    echo "Unrecognized EXE directory!"
    exit
fi

echo "INP="${INP}
echo "OUT="${OUT}
echo ""

echo "CMD=${CMD}"
echo ""

eval $CMD

echo "done merge_builds..."

# End of file
