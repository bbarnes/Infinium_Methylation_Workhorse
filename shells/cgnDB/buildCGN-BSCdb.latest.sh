

# Original build as of 20022021 with Top Sequence Output::
#  RUN=GRCh36-38.mm10.FWD_SEQ

# New build with Forward Sequences::
# RUN=GRCh36-GRCh37-GRCh38-GRCm10.12032021
# RUN=GRCh36-GRCh37-GRCh38-GRCm10.13032021
RUN=GRCh36-GRCh37-GRCh38-GRCm10.14032021_test

TOP=/illumina/scratch/darkmatter/Projects/dbCGN
EXE=${TOP}/scripts/cgnDB/buildCGN-BSCdb.pl
OUT=${TOP}/data/cgnDB/${RUN}/design-input

#
# Need an option for reading in csv version...
#
# CGN=${TOP}/manifests/manifest.cgn.hash
# CGN=/illumina/scratch/darkmatter/dbCGN/manifests/manifest.cgn.hash
# CGN=${TOP}/data/cgnDB/GRCh36-38.mm10/pre-assigned.cgnTop.hash.csv.gz
#
CGN="${TOP}/data/cgnDB/GRCh36-38.mm10/canonical-12032021.cgn-top-src.hash.csv.gz"

GEN38=/illumina/scratch/darkmatter/iGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/GRCh38.genome.fa.gz
GEN37=/illumina/scratch/darkmatter/iGenomes/Homo_sapiens/NCBI/GRCh37/Sequence/WholeGenomeFasta/GRCh37.genome.fa.gz
GEN36=/illumina/scratch/darkmatter/iGenomes/Homo_sapiens/NCBI/GRCh36/Sequence/WholeGenomeFasta/GRCh36.genome.fa.gz
GEN10=/illumina/scratch/darkmatter/iGenomes/Mus_musculus/NCBI/GRCm10/Sequence/WholeGenomeFasta/GRCm10.genome.fa.gz

mkdir -p ${OUT}

${EXE} -v 3 -out ${OUT} \
    -cgn ${CGN} \
    -mod a -name ${RUN} \
    -src GRCh38 -f ${GEN38}

#    -src GRCh37 -f ${GEN37} \
#    -src GRCh36 -f ${GEN36} \
#    -src GRCm10 -f ${GEN10}

echo "done."

exit




DES38=${TOP}/cgn/GRCh36-38.mm10/GRCh38.improbeDesignInput.tsv.gz
DES37=${TOP}/cgn/GRCh36-38.mm10/GRCh37.improbeDesignInput.tsv.gz
DES36=${TOP}/cgn/GRCh36-38.mm10/GRCh36.improbeDesignInput.tsv.gz
DES10=${TOP}/cgn/GRCh36-38.mm10/mm10.improbeDesignInput.tsv.gz

NAM=GRCh37-38
OUT=${TOP}/cgn/${NAM}
${EXE} -v 3 -out ${OUT} \
    -mod w -name ${NAM} \
    -d ${DES37} \
    -d ${DES38} \

NAM=GRCh38
OUT=${TOP}/cgn/${NAM}
${EXE} -v 3 -out ${OUT} \
    -mod w -name ${NAM} \
    -d ${DES38}

NAM=GRCh37
OUT=${TOP}/cgn/${NAM}
${EXE} -v 3 -out ${OUT} \
    -mod w -name ${NAM} \
    -d ${DES37}

NAM=GRCh36
OUT=${TOP}/cgn/${NAM}
${EXE} -v 3 -out ${OUT} \
    -mod w -name ${NAM} \
    -d ${DES36}

NAM=mm10
OUT=${TOP}/cgn/${NAM}
${EXE} -v 3 -out ${OUT} \
    -mod w -name ${NAM} \
    -d ${DES10}

## end of file
