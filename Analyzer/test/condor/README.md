# Condor production on uaf

```
cmsenv; scram b enable-multi-targets; scram b -j8   # the worker pool mixes AVX2 and pre-AVX2 CPUs
mkdir -p run && cd run && mkdir -p logs output
../make_ship.sh
../make_jobs_data.sh Muon0 20 "saveTrigNames=False;tpOnly=True" > jobs.txt
cp ../run_mcpana.sh ../mcpana.sub .
condor_submit mcpana.sub joblist=jobs.txt
```
Each job line is `<name> <maxEvents> <cfg args separated by ;> <comma-separated /store files> [cfg]`
(`cfg` defaults to `mcpAnalyzer_cfg.py`; use `mcpZProbe_cfg.py` for RECO Z probes). Join several input files with `+`, not `,`: condor splits queue columns on commas;
outputs land in `output/<name>.root` with the cmsRun log `output/cmsrun_<name>.log` (last line `exit <code>`).
Workers need the rhel8 singularity image (set in `mcpana.sub`), and read inputs over xrootd.
