#!/bin/sh

RSCRIPT_A=/usr/bin/Rscript
RSCRIPT_B=/usr/local/bin/Rscript

if [ -e ${RSCRIPT_A} ]; then

    RSCRIPT=${RSCRIPT_A}

elif [ -e ${TOP_LIX} ]; then

    RSCRIPT=${RSCRIPT_B}

else
    echo "Unrecognized RSCRIPT directory!"
    exit
fi

${RSCRIPT} \
    /repo/Infinium_Methylation_Workhorse/scripts/R/swifthoof/swifthoof_main.R \
    --Rscript /usr/local/bin/Rscript \
    -i . -o /output \
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
    --parallel \
    --verbose=30


#
# --datDir /repo/workhorse/dat \
# --platform=EPIC \
# --manifest=B4 \

# End of file
