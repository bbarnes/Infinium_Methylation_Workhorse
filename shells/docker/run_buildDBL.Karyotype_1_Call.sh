#!/bin/sh

/usr/bin/Rscript \
    /repo/Infinium_Methylation_Workhorse/scripts/R/analysis/build_models.R \
    --Rscript /usr/bin/Rscript \
    -m . -o /output \
    --buildDbl \
    --runName=r1 \
    --classVar=Karyotype_1_Call \
    --trainClass=XaXi,XaY \
    --samplePvalName=Poob_Pass_0_Perc \
    --samplePvalPerc=96 \
    --alphaMin=0.0 \
    --alphaMax=1.0 \
    --alphaInc=0.5 \
    --seeds=13,42 \
    --featuresDml=100,1000 \
    --lociBetaKey=ind_beta \
    --lociPvalKey=i_poob \
    --lociPvalMin=0.9 \
    --parallel \
    --cluster \
    --verbose=6


#    --buildDbl \
#    --buildModels \


# End of file
