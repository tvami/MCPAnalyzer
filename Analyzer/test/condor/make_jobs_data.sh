#!/bin/bash
# job list: one file per job, every Nth file of each 2024 MINIv6NANOv15 dataset of a PD
# usage: make_jobs_data.sh <PD, e.g. Muon0> <N> [extra cfg args, ';'-separated] > jobs.txt
PD=$1; N=$2; EXTRA=${3:-saveTrigNames=False}
for ds in $(dasgoclient -query="dataset=/$PD/Run2024*MINIv6NANOv15*/MINIAOD"); do
  era=$(echo $ds | cut -d/ -f3 | cut -c8); proc=$(echo $ds | cut -d/ -f3 | cut -d- -f2- | tr -d '-')
  dasgoclient -query="file dataset=$ds" | sort | awk -v e=$era -v p=$proc -v pd=$PD -v n=$N -v x="$EXTRA" \
    'NR%n==1 || n==1 {printf "%s_%s_%s_%05d -1 isData=True;lumiMask=json/2024%s_Golden.json;%s %s\n", pd, e, p, NR, e, x, $1}'
done
