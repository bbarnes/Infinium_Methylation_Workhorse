#!/bin/bash

PROGRAM="Infinium_Methylation_Workhorse/scripts/R/workhorse/workhorse_main.R"

DOC_VER="v.1.55"
DOC_NAME="Infinium_Methylation_Workhorse_Centos"
DOC_IMAGE="bbarnesimdocker/im_workhorse:${DOC_NAME}.${DOC_VER}"
DOC_SHELL="run_cmd.sh"

VERBOSE="5"
SPECIES="Homo_sapiens"

TOP_DIR_A="/illumina/scratch/darkmatter"
TOP_DIR_B="/Users/bretbarnes/Documents"

SYSTEM="MAC"

ORD_CSV="McMaster_CpG_DesignFile_v4.csv.gz"
MAT_TSV="20532820_probes1.match.gz,20532820_probes2.match.gz"
AQP_TSV="20051339_A_ProductQC.txt.gz"

RELOAD="--reload"
RELOAD=""

PARALLEL=""
PARALLEL="--parallel"

if [ -e ${TOP_DIR_A} ]; then
    TOP_DIR=${TOP_DIR_A}

    SIG_IMAGE="/illumina/scratch/darkmatter/docker/software/${DOC_IMAGE}.sif"

    BSP_DIR="/illumina/scratch/methylation/programs/..."
    BSP_EXE="bsmap"

    # MAN_LDIR="/illumina/scratch/darkmatter/tools/Infinium_Methylation_Workhorse/dat/manifest/core"
    # MAN_SDIR="-B ${MAN_LDIR:/tmp}"

elif [ -e ${TOP_DIR_B} ]; then
    TOP_DIR=${TOP_DIR_B}

    BSP_DIR="${TOP_DIR}/tools/programs/BSMAPz"
    BSP_EXE="bsmap"

else
    echo "Unrecognized Rscript EXE!"
    exit
fi

RUN_NAME="McMaster10Kselection"

TAG_MAP_TSV="GRCh37.chr-pos-srd.slim.cgn-sorted.txt.gz"
BSP_MAP_TSV="GRCh37.chr-pos-srd.slim.pos-sorted.txt.gz"

OPT_STR="${PROGRAM} \
  --run_name=${RUN_NAME} \
  --platform=MCM \
  --version=v1 \
  --genome_build=GRCh37 \
  --ord_csv=${ORD_CSV} \
  --mat_tsv=${MAT_TSV} \
  --aqp_tsv=${AQP_TSV} \
  --canonical_cgn_csv=canonical.cgn-top-grp.csv.gz \
  --tag_map_tsv=${TAG_MAP_TSV} \
  --bsp_map_tsv=${BSP_MAP_TSV} \
  --trackTime \
  --Rscript=Rscript \
  --bsmap_exe=${BSP_EXE} \
${RELOAD} ${PARALLEL} --verbose=${VERBOSE} "

# --memory-swap="[memory_limit]"
# --memory=${MEM}
MEM="16g"

ORD_DIR="${TOP_DIR}/data/CustomContent/McMaster/McMaster10Kselection/AQP.v2/order"
MAT_DIR="${TOP_DIR}/data/CustomContent/McMaster/McMaster10Kselection/AQP.v2/match"
AQP_DIR="${TOP_DIR}/data/CustomContent/McMaster/McMaster10Kselection/AQP.v2/aqp"

OUT_DIR="${TOP_DIR}/scratch/docker-${DOC_VER}"

GEN_DIR="${TOP_DIR}/data/imGenomes/${SPECIES}/NCBI"
MAN_DIR="${TOP_DIR}/data/manifests"
IMP_DIR="${TOP_DIR}/data/improbe"
ANN_DIR="${TOP_DIR}/data/annotation"

SEQ_DIR="${IMP_DIR}/scratch/cgnDB/dbSNP_Core4/design-output/prbs-p49-split"
POS_DIR="${IMP_DIR}/scratch/cgnDB/dbSNP_Core4/design-input/min"
CAN_DIR="${IMP_DIR}/data/improbe/scratch/cgnDB/dbSNP_Core4/design-input"

# If you just want to verify docker works::
#  print the help command::
#
# OPT_STR="${PROGRAM} --help "

mkdir -p ${OUT_DIR}

echo "Local Mounted Directoreis::"
echo "  OUT_DIR = ${OUT_DIR}"
echo "  BSP_DIR = ${BSP_DIR}"
echo
echo "  ORD_DIR = ${ORD_DIR}"
echo "  MAT_DIR = ${MAT_DIR}"
echo "  AQP_DIR = ${AQP_DIR}"
echo
echo "  GEN_DIR = ${GEN_DIR}"
echo "  MAN_DIR = ${MAN_DIR}"
echo "  IMP_DIR = ${IMP_DIR}"
echo "  ANN_DIR = ${ANN_DIR}"
echo
echo "  SEQ_DIR = ${SEQ_DIR}"
echo "  POS_DIR = ${POS_DIR}"
echo "  CAN_DIR = ${CAN_DIR}"
echo

echo "Input file names (order, match, aqp/pqc)::"
echo "  ORD_CSV = ${ORD_CSV}"
echo "  MAT_TSV = ${MAT_TSV}"
echo "  AQP_TSV = ${AQP_TSV}"
echo

if [ -e ${TOP_DIR_A} ]; then

    SCMD="singularity exec \
	 -B ${ORD_DIR}:/order \
	 -B ${BSP_DIR}:/bsp \
	 -B ${MAT_DIR}:/match \
	 -B ${AQP_DIR}:/aqp \
	 -B ${OUT_DIR}:/output \
	 -B ${SEQ_DIR}:/sequence \
	 -B ${POS_DIR}:/coordinate \
	 -B ${CAN_DIR}:/canonical \
	 -B ${GEN_DIR}:/genome \
	 -B ${IMP_DIR}:/improbe \
	 -B ${ANN_DIR}:/annotation \
	 -B ${MAN_DIR}:/manifest \
      	 -B ${SCRATCH_ROOT}:${SCRATCH_ROOT} \
      	 --workdir ${SCRATCH_ROOT} \
      	 ${SIG_IMAGE} ${DOC_SHELL} ${OPT_STR}"

elif [ -e ${TOP_DIR_B} ]; then

    RUN_CMD="docker run -i --rm -w /work \
             -v ${ORD_DIR}:/order \
     	     -v ${MAT_DIR}:/match \
     	     -v ${AQP_DIR}:/aqp \
     	     -v ${OUT_DIR}:/output \
     	     -v ${SEQ_DIR}:/sequence \
     	     -v ${POS_DIR}:/coordinate \
     	     -v ${CAN_DIR}:/canonical \
     	     -v ${GEN_DIR}:/genome \
     	     -v ${IMP_DIR}:/improbe \
     	     -v ${ANN_DIR}:/annotation \
     	     -v ${MAN_DIR}:/manifest \
     	     ${DOC_IMAGE} ${DOC_SHELL} ${OPT_STR}"

else
    echo "Unrecognized Rscript EXE!"
    exit
fi

echo
echo "RUN_CMD = ${RUN_CMD}"
echo

${RUN_CMD}

echo "Done: $0"

# End of file
