#!/bin/bash
NREAS=10

RATE=$(basename "`pwd`")
POL=$(pwd | xargs dirname | xargs basename)
BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

for i in $(seq -f %02g 1 $NREAS); do
    mkdir -p $i
    cd $i
    pwd
    relpath=$(realpath --relative-to="." "$BASE")
    sed -e "s/--job-name=*/--job-name=$POL\_$RATE\_rea\_$i/g" \
        -e "s|BASE|$relpath|g" $BASE/srun_pull_template.sh > srun.sh
    sbatch srun.sh
    cd ..
done
