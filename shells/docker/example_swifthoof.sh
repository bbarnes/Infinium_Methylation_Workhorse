
VER=v.3.32

EXP=COVIC-Set1-15052020
EXP=CNTL-Samples_VendA_10092020
DAT=/Users/bretbarnes/Documents/data/idats/idats_${EXP}
OUT=/Users/bretbarnes/Documents/scratch/docker-${VER}/swifthoof/${EXP}

SRC=bbarnesimdocker/im_workhorse:Infinium_Methylation_Workhorse.${VER}
CMD=run_swifthoof.sh

docker run -i --rm -v ${DAT}:/input -v ${OUT}:/output -w /work ${SRC} ${CMD} \
       --workflows=ind  \
       --minNegPval=0.02 --minOobPval=0.1 --minNegPerc=98 --minOobPerc=90 --minDeltaBeta=0.2 \
       --verbose=8 \
       --write_call \
       --single

echo "swifthoof done..."

# End of file
