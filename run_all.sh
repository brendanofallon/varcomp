#!/bin/bash


VERSION="v8"
ZYG="het"
RUNDIR="rundirs_${ZYG}_${VERSION}"
MAX_JOBS=3

CMD="python /slc-ngs/projects/varcomp/vcomp/injectvar.py -c ../comp.conf -v ../clinvar_{1}.vcf --het > results_{1}_${ZYG}_${VERSION}.txt"



mkdir $RUNDIR
cd  $RUNDIR
parallel --jobs=$MAX_JOBS $CMD ::: {1..7}


