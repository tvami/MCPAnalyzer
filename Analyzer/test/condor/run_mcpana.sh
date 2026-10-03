#!/bin/bash
# args: outname maxEvents extra_cfg_args inputFile(s, comma separated) [cfg, default mcpAnalyzer_cfg.py]
mkdir -p output
OUT=$1; NEV=$2; EXTRA=$3; IN=$4; CFGNAME=${5:-mcpAnalyzer_cfg.py}
trap 'cp -f *.log output/ 2>/dev/null' EXIT
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el8_amd64_gcc12
scram project CMSSW CMSSW_15_0_19_patch2 > /dev/null || exit 3
tar xzf mcpana_ship.tgz -C CMSSW_15_0_19_patch2
tar xzf golden24.tgz 2>/dev/null
cd CMSSW_15_0_19_patch2/src && eval `scramv1 runtime -sh` && cd ../..
FILES=$(echo $IN | tr ',+' '\n\n' | sed 's|^/store|root://cms-xrd-global.cern.ch//store|' | paste -sd,)
cmsRun $CMSSW_BASE/src/MCPAnalyzer/Analyzer/test/$CFGNAME inputFiles=$FILES maxEvents=$NEV outputFile=output/$OUT.root $(echo $EXTRA | tr ';' ' ') > cmsrun_$OUT.log 2>&1
rc=$?; echo "exit $rc" >> cmsrun_$OUT.log; exit $rc
