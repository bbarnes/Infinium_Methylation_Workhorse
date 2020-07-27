#!/bin/sh

/usr/bin/Rscript \
    /repo/Infinium_Methylation_Workhorse/scripts/R/analysis/merge_builds.R \
    --Rscript /usr/bin/Rscript \
    -b . -o /output \
    --runName=r1 \
    --platform=EPIC \
    --version=B4 \
    --classVar=Karyotype_1_Call \
    --parallel \
    --cluster \
    --verbose=6

# End of file
