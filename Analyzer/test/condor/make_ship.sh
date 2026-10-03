#!/bin/bash
# pack the local build (incl. the x86-64-v2 variant) and the 2024 Golden JSONs for condor
# needs: cmsenv in CMSSW_15_0_19_patch2, `scram b enable-multi-targets` done once
set -e
tar czf mcpana_ship.tgz -C $CMSSW_BASE lib/$SCRAM_ARCH cfipython/$SCRAM_ARCH/MCPAnalyzer src/MCPAnalyzer/Analyzer/test/mcpZProbe_cfg.py \
    src/MCPAnalyzer/Analyzer/test/mcpAnalyzer_cfg.py
mkdir -p json && cp /cvmfs/cms-griddata.cern.ch/cat/metadata/DC/Collisions24/latest/2024{C,D,E,F,G,H,I}_Golden.json json/
tar czf golden24.tgz json && rm -rf json
