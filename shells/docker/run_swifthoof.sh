#!/bin/sh

EXE_NAME=Infinium_Methylation_Workhorse/scripts/R/swifthoof/swifthoof_main.R
EXE_A=/repo/${EXE_NAME}
EXE_B=/Users/bretbarnes/Documents/tools/${EXE_NAME}

RSCRIPT_A=/usr/local/bin/Rscript
RSCRIPT_B=/usr/local/bin/Rscript

if [ -e ${EXE_A} ]; then
    EXE=${EXE_A}

    INP=/input
    OUT=/output

    echo "Docker Run..."
    
elif [ -e ${EXE_B} ]; then
    EXE=${EXE_B}

    echo "Local Run..."
    
else
    echo "Unrecognized EXE directory!"
    exit
fi

echo "RSC="${RSCRIPT}
echo "EXE="${EXE}
echo "INP="${INP}
echo "OUT="${OUT}
echo ""

echo ${RSCRIPT} ${EXE} --Rscript ${RSCRIPT} -i ${INP} -o ${OUT} $@

${RSCRIPT} ${EXE} --Rscript ${RSCRIPT} --idatsDir ${INP} --outDir ${OUT} $@

echo "done."

exit

${RSCRIPT} ${EXE} \
    --Rscript ${RSCRIPT} \
    -i ${INP} -o ${OUT} \
    --workflows="i,ind" \
    --writeCalls \
    --writeSsheet


    
#    --minNegPval=0.02 \
#    --minOobPval=0.1 \
#    --minNegPerc=98 \
#    --minOobPerc=90 \
#    --minDeltaBeta=0.2 \
#    --percisionBeta=4 \
#    --percisionPval=6 \
#    --verbose=${VER}

# End of file
