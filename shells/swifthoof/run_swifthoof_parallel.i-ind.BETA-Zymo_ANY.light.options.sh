#!/bin/bash

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 outTopDir idatsDir experimentName [options]"
    exit 1
fi
outTop=$1
idatsDir=$2
EXP_NAME=$3
optSTR=$4

# optSTR=echo ${optSTR} | sed 's/([^,]+)/--$1 /gi'
# optSTR2=${optSTR//,+$/ /}
optSTR2=${optSTR//,/ }

# outSuffix=$3
# EXP_NAME=`basename $idatsDir`
# EXP_NAME=$(sed 's/^idats_//' <<< "$EXP_NAME")

verbose=3

# Platform/Method Options::
# platform="HM450"
# platform="EPIC"

# manifest="B2"
# manifest="B4"
# manifest="C0"

# Run Options::
fresh=true
buildSubDir=true
# autoDetect=true
autoDetect=false
workflows='i,ind'

# Output Options::
loadIdat=true
saveIdat=false

skipSwap=false
loadSsets=false
saveSsets=false
saveRawSset=false

writeSset=false
writeCalls=true
writeSsheet=true
writeAuto=false

# Optional Files::
subManifest=false

# Threshold Options::
minNegPval=0.02
minOobPval=0.1
minNegPerc=98
minOobPerc=90
minDeltaBeta=0.2

percisionBeta=4
percisionPval=6

# Parallel/Cluster Options::
# single=true
# parallel=true
# cluster=true

single=false
parallel=false
cluster=false

# Plotting Options::
plotSset=false
plotCalls=false
plotAuto=false

plotFormat="pdf"
plotFormat="png"

dpi=72
dpi=120
plotMax=10000
plotSub=5000

# Standard Options:: DON'T CHANGE::
#
prgmTop="Infinium_Methylation_Workhorse"
prgmDir="swifthoof"
prgmTag="swifthoof_main"

TOP_MAC1=/Users/bbarnes/Documents/Projects/methylation/tools
TOP_MAC2=/Users/bretbarnes/Documents/tools
TOP_LIX1=/illumina/scratch/darkmatter/tools

if [ -e ${TOP_MAC1} ]; then
    TOP=${TOP_MAC1}
    RSCRIPT=/usr/local/bin/Rscript

elif [ -e ${TOP_MAC2} ]; then
    TOP=${TOP_MAC2}
    RSCRIPT=/usr/local/bin/Rscript

elif [ -e ${TOP_LIX1} ]; then
    TOP=${TOP_LIX1}

    CONDA=conda_4.6.8
    # CONDA=Anaconda2-2019.10-Linux-x86_64
    # CONDA=Anaconda3-2019.10-Linux-x86_64
    # CONDA=Anaconda3-2019.10-Linux-x86_64.2020
    RSCRIPT=/illumina/scratch/darkmatter/thirdparty/${CONDA}/bin/Rscript

else
    echo "Unrecognized top directory!"
    exit
fi

SRC=${TOP}/${prgmTop}
DAT=${SRC}/dat
EXE=${SRC}/scripts/R/${prgmDir}/${prgmTag}.R

CMD=${RSCRIPT}" "${EXE}
CMD+=" --"Rscript=${RSCRIPT}

# Directories::
#  outDir=${TOP}/builds/${prgmTag}${outSuffix}/${EXP_NAME}
outDir=${outTop}/${prgmTag}${outSuffix}/${EXP_NAME}

auto_sam_csv=${DAT}/ref/AutoSampleDetection_EPIC-B4_8x1_pneg98_Median_beta_noPval_BETA-Zymo.csv.gz
auto_sam_csv=${DAT}/ref/AutoSampleDetection_EPIC-B4_8x1_pneg98_Median_beta_noPval_BETA-Zymo_Mean-COVIC-90-NP-ind_negs-0.02.csv.gz
auto_sam_csv=${DAT}/ref/AutoSampleDetection_EPIC-B4_8x1_pneg98_Median_beta_noPval_BETA-Zymo_Mean-COVIC-280-NP-ind_negs-0.02.csv.gz

# idatsDir=NULL
CMD+=" --auto_sam_csv"=${auto_sam_csv}
CMD+=" --outDir"=${outDir}
CMD+=" --idatsDir"=${idatsDir}

# Optional Files::
if [ "${subManifest}" = true ]; then
    CMD+=" --subManifest"
fi
# CMD+=" --manifestPath"=${manifestPath}
# CMD+=" --addressPath"=${addressPath}

# Platform/Method Options::
# CMD+=" --platform"=${platform}
# CMD+=" --manifest"=${manifest}

# Run Options::
if [ "${fresh}" = true ]; then
    CMD+=" --fresh"
fi
if [ "${buildSubDir}" = true ]; then
    CMD+=" --buildSubDir"
fi
if [ "${autoDetect}" = true ]; then
    CMD+=" --autoDetect"
fi
CMD+=" --workflows"=${workflows}

# Output Options::
if [ "${loadIdat}" = true ]; then
    CMD+=" --loadIdat"
fi
if [ "${saveIdat}" = true ]; then
    CMD+=" --saveIdat"
fi

if [ "${skipSwap}" = true ]; then
    CMD+=" --skipSwap"
fi
if [ "${loadSsets}" = true ]; then
    CMD+=" --loadSsets"
fi
if [ "${saveSsets}" = true ]; then
    CMD+=" --saveSsets"
fi
if [ "${saveRawSset}" = true ]; then
    CMD+=" --saveRawSset"
fi

if [ "${writeSset}" = true ]; then
    CMD+=" --writeSset"
fi
if [ "${writeCalls}" = true ]; then
    CMD+=" --writeCalls"
fi
if [ "${writeSsheet}" = true ]; then
    CMD+=" --writeSsheet"
fi
if [ "${writeAuto}" = true ]; then
    CMD+=" --writeAuto"
fi

# Threshold Options::
CMD+=" --minNegPval"=${minNegPval}
CMD+=" --minOobPval"=${minOobPval}
CMD+=" --minNegPerc"=${minNegPerc}
CMD+=" --minOobPerc"=${minOobPerc}
CMD+=" --minDeltaBeta"=${minDeltaBeta}
CMD+=" --percisionBeta"=${percisionBeta}
CMD+=" --percisionPval"=${percisionPval}

# Parallel/Cluster Options::
if [ "${single}" = true ]; then
    CMD+=" --single"
fi
if [ "${parallel}" = true ]; then
    CMD+=" --parallel"
fi
if [ "${cluster}" = true ]; then
    CMD+=" --cluster"
fi

# Plotting Options::
if [ "${plotSset}" = true ]; then
    CMD+=" --plotSset"
fi
if [ "${plotCalls}" = true ]; then
    CMD+=" --plotCalls"
fi
if [ "${plotAuto}" = true ]; then
    CMD+=" --plotAuto"
fi

CMD+=" --plotFormat"=${plotFormat}
CMD+=" --dpi"=${dpi}
CMD+=" --plotMax"=${plotMax}
CMD+=" --plotSub"=${plotSub}

# Verbosity Options::
CMD+=" --"verbose=${verbose}
CMD+=" "${optSTR2}

# mkdir -p ${outDir}
echo ${CMD}

${CMD}

## End of file