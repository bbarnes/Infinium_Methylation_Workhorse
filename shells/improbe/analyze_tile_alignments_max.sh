#!/bin/bash

if [ "$#" -lt 7 ]; then
    echo "Usage: $0 aln_dir outDir runName platform version build genome MAX"
    exit 1
fi
aln_dir=$1
outDir=$2
runName=$3

platform=$4
version=$5
build=$6
genome=$7
max=$8

verbose=3

# Parallel/Cluster Options::
# single=true
single=false
parallel=true
cluster=true

# Program Variables::
#
prgmTop="Infinium_Methylation_Workhorse"
prgmDir="probe_design"
prgmTag="analyze_tile_alignments"

TOP_MAC=/Users/bbarnes/Documents/Projects/methylation/tools
TOP_LIX=/illumina/scratch/darkmatter/Projects/COVIC/tools

if [ -e ${TOP_MAC} ]; then
    TOP=${TOP_MAC}
    CONDA=mac
    RSCRIPT=/usr/local/bin/Rscript
elif [ -e ${TOP_LIX} ]; then
    TOP=${TOP_LIX}
    CONDA=conda_4.6.8
    # CONDA=Anaconda2-2019.10-Linux-x86_64
    # CONDA=Anaconda3-2019.10-Linux-x86_64
    RSCRIPT=/illumina/scratch/darkmatter/thirdparty/${CONDA}/bin/Rscript
else
    echo "Unrecognized top directory!"
    exit
fi

SRC=${TOP}/${prgmTop}
EXE=${SRC}/scripts/R/${prgmDir}/${prgmTag}.R

CMD=${RSCRIPT}" "${EXE}
CMD+=" --"Rscript=${RSCRIPT}

CMD+=" --outDir"=${outDir}
CMD+=" --runName"=${runName}
CMD+=" --aln_dir"=${aln_dir}
CMD+=" --platform"=${platform}
CMD+=" --version"=${version}
CMD+=" --build"=${build}
CMD+=" --genome"=${genome}

CMD+=" --max"=${max}

if [ "${single}" = true ]; then
    CMD+=" --single"
fi
if [ "${parallel}" = true ]; then
    CMD+=" --parallel"
fi
if [ "${cluster}" = true ]; then
    CMD+=" --cluster"
fi


# Verbosity Options::
CMD+=" --"verbose=${verbose}

# mkdir -p ${outDir}
echo ${CMD}

${CMD}

## End of file
