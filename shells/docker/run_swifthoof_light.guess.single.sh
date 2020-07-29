#!/bin/sh

/usr/bin/Rscript \
    /repo/Infinium_Methylation_Workhorse/scripts/R/swifthoof/swifthoof_main.R \
    --Rscript /usr/bin/Rscript \
    -i . -o /output \
    --workflows=ind \
    --lightFootPrint \
    --writeCalls \
    --writeSsheet \
    --sigs_sum_field=avg \
    --minNegPval=0.02 \
    --minOobPval=0.1 \
    --minNegPerc=98 \
    --minOobPerc=90 \
    --minDeltaBeta=0.2 \
    --percisionBeta=4 \
    --percisionPval=6 \
    --verbose=3 \
    --single


#    --parallel \
# --datDir /repo/workhorse/dat \
# --platform=EPIC \
# --manifest=B4 \

# End of file