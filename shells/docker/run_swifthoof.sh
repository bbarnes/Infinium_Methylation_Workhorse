#!/bin/sh

EXE_A=/repo/Infinium_Methylation_Workhorse/scripts/R/swifthoof/swifthoof_main.R
EXE_B=/Users/bretbarnes/Documents/tools/Infinium_Methylation_Workhorse/scripts/R/swifthoof/swifthoof_main.R

RSCRIPT_A=/usr/bin/Rscript
RSCRIPT_B=/usr/local/bin/Rscript

if [ -e ${EXE_A} ]; then
    EXE=${EXE_A}

    INP=/input
    OUT=/output
    
elif [ -e ${EXE_B} ]; then
    EXE=${EXE_B}
    
else
    echo "Unrecognized EXE directory!"
    exit
fi

if [ -e ${RSCRIPT_A} ]; then
    RSCRIPT=${RSCRIPT_A}

elif [ -e ${RSCRIPT_B} ]; then
    RSCRIPT=${RSCRIPT_B}

    if [ "$#" -lt 2 ]; then
	echo "Usage: $0 input output [options]"
	exit 1
    fi

    INP=$1
    OUT=$2

else
    echo "Unrecognized RSCRIPT directory!"
    exit
fi

echo "RSC="${RSCRIPT}
echo "EXE="${EXE}
echo "INP="${INP}
echo "OUT="${OUT}
echo ""

echo ${RSCRIPT} ${EXE} -Rscript ${RSCRIPT} -i ${INP} -o ${OUT} $@

${RSCRIPT} ${EXE} -Rscript ${RSCRIPT} -i ${INP} -o ${OUT} $@

echo "done."

exit

${RSCRIPT} ${EXE} \
    --Rscript /usr/local/bin/Rscript \
    -i ${INP} -o ${OUT} \
    --workflows="i,ind" \
    --writeCalls \
    --writeSsheet \
    --minNegPval=0.02 \
    --minOobPval=0.1 \
    --minNegPerc=98 \
    --minOobPerc=90 \
    --minDeltaBeta=0.2 \
    --percisionBeta=4 \
    --percisionPval=6 \
    --verbose=${VER}

# End of file
