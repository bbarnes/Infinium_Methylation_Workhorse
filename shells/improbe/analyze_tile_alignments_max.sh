
MAX=20
MAX=3620
MAX=100

TOP=/illumina/scratch/darkmatter/Projects/COVIC
ALN=/illumina/scratch/darkmatter/Projects/COVIC/scratch/index/tile_main/EPIC/SARS-CoV-2/MN908947/COVIC/align

# TOP=/Users/bbarnes/Documents/Projects/methylation
# ALN=/Users/bbarnes/Documents/Projects/methylation/scratch/small.index/tile_main/EPIC/SARS-CoV-2/MN908947/COVIC/align

EXE=${TOP}/tools/Infinium_Methylation_Workhorse/shells/improbe/analyze_tile_alignments_max.sh
FAS=${TOP}/fas/nCoV_Wuhan_Sequence_MN908947.3.fa.gz
FAS=${TOP}/fas/nCoV_Wuhan_Sequence_MN908947.3.fa
OUT=${TOP}/scratch

# GEN=/illumina/scratch/darkmatter/bscGenomes/iGenomes/hg38/hg38.genome.fa,/illumina/scratch/darkmatter/bscGenomes/iGenomes/hg37/hg37.genome.fa,/illumina/scratch/darkmatter/Projects/COVIC/fas/nCoV_Wuhan_Sequence_MN908947.3.fa,/illumina/scratch/darkmatter/Projects/COVIC/fas/ncbi-sequences.fas

# HUM_GEN=/illumina/scratch/darkmatter/bscGenomes/iGenomes/hg38/hg38.genome.fa,/illumina/scratch/darkmatter/bscGenomes/iGenomes/hg37/hg37.genome.fa
# GEN=/illumina/scratch/darkmatter/Projects/COVIC/fas/nCoV_Wuhan_Sequence_MN908947.3.fa,/illumina/scratch/darkmatter/Projects/COVIC/fas/ncbi-sequences.fas

# Index Test::
#
GEN=/illumina/scratch/darkmatter/Projects/COVIC/fas/ncbi_index

OUT=${TOP}/scratch

RNAME="COVIC"
PLAT="EPIC"
VER="SARS-CoV-2"
BUILD="MN908947"

echo ${EXE} ${ALN} ${OUT} ${RNAME} ${PLAT} ${VER} ${BUILD} ${GEN} ${MAX}

# 0       1      2      3        4       5       6       7     8
${EXE} ${ALN} ${OUT} ${RNAME} ${PLAT} ${VER} ${BUILD} ${GEN} ${MAX}

# End of file
