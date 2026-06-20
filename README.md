# MCPAnalyzer

Self-contained package to characterize **HSCP multi-charged particles (MCP)** on MiniAOD and to
develop pixel/strip discriminators against SM background. The MCP is `|pdgId| = 10000200`.

Runs inside `cmsRun` (needs a 2024 MC GlobalTag + geometry + the pixel CPE). It reads the
`pat::IsolatedTrack` collection and the `DeDxHitInfo` association attached to it — no
`HSCParticleProducer` needed — and writes two flat TTrees (`tracks`, `gen`) via `TFileService`.

## Build & run
```
cd CMSSW_15_0_19_patch2/src && cmsenv && scram b -j8
cmsRun MCPAnalyzer/Analyzer/test/mcpAnalyzer_cfg.py \
    inputFiles=file:<MiniAOD.root> maxEvents=-1 outputFile=out.root
```
Key config parameters (`test/mcpAnalyzer_cfg.py`): `gtag` (default `150X_mcRun3_2024_realistic_v2`),
`pixelCPE` (`PixelCPETemplateReco`).

On a background sample (e.g. QCD) there is no MCP, so every track has `genMatched=0`; the
discriminant branches are still filled for all tracks.

### Producing the ntuples used in the studies
All runs are from `CMSSW_15_0_19_patch2/src` after `cmsenv`. Output goes to `~/HSCP-MCP/`.

**Signal** (M=2000 charge ladder; the MiniAOD lives next to the GEN-SIM under the CRAB tree):
```
B=/ceph/cms/store/user/tvami/HSCP_MCP
for Q in 6 30 90 ; do          # Q-tag = 3 x charge(e); 6->2e, 30->10e, 90->30e
  f=$(ls $B/HSCP-MCP_M-2000-Q-${Q}_GEN-SIM_2024_v2_13p6TeV/crab_HSCP-MCP_M-2000-Q-${Q}_MiniAOD_2024_v2/*/0000/*.root | head -1)
  cmsRun MCPAnalyzer/Analyzer/test/mcpAnalyzer_cfg.py \
      inputFiles=file:$f maxEvents=-1 outputFile=$HOME/HSCP-MCP/mcp_M2000_Q${Q}_v2.root
done
```

**Background** (muon-enriched QCD, `RunIII2024Summer24MiniAODv6-150X`, one file copied into `src/`):
```
cmsRun MCPAnalyzer/Analyzer/test/mcpAnalyzer_cfg.py \
    inputFiles=file:000b4dc9-1495-498c-b742-4ca912780503.root \
    maxEvents=-1 outputFile=$HOME/HSCP-MCP/mcp_QCD_Pt50to80_v2.root
```
The QCD ntuple is the per-track discriminant-shape reference (its tracks all have `genMatched=0`).

## Discriminating variables (in the `tracks` tree)

### Charge probability — `probQ_pixel`
- **`probQ_pixel`** — our port of the standard pixel `probQ`: per pixel hit, the SiPixel template
  (via the CPE) gives `probabilityQ`, the p-value of the cluster charge against the expected MIP
  charge for the track's angle; the per-hit values are combined with the Poisson series
  `combineProbs`. The companion **`probXY_pixel`** is the gated cluster-shape probability
  (`probQ<0.8 && probXY>0`), exactly as in the standard tool. `probQ_pixelNoL1` excludes BPIX L1.
  - **Never combine a per-hit `probQ == 0`.** A per-hit value at/below `kProbHitFloor = 1e-6` is a
    saturated / off-template cluster (a *missing* measurement, not a "probability 0" one) and is
    skipped from the product. This stops a single pathological hit from underflowing the float
    product to a fake track-level `0`. If **no** usable charge hit survives, the track value is set
    to **`-1` (undefined)**, not `0` — otherwise the most MIP-like tracks (every hit `probQ>0.8`,
    so `numRecHits==0` → old `combineProbs(1,0)==0`) would pile up at the signal-like cap. Plotting
    code must drop `probQ < 0` (undefined) and may keep an exact `0` only as genuine underflow.

### Energy loss — `ih_pixel`, `ih_strip`, `ih_full`
Harmonic-mean dE/dx (MeV/cm), pixel-only / strip-only / combined, with `nom_*` measurement counts.
NB: `probQ_pixel` and `ih_pixel` are **not independent** — both are the pixel cluster charge.

### Pixel cluster morphology — `pixClSizeX` and the angle-corrected `pixSizeXresidual`
Raw cluster sizes (`pixClSize`, `pixClSizeX`, `pixClSizeY`, `pixClSizeMax`, `pixClCharge`) are stored,
but the useful one is the **angle-corrected r-φ cluster size**.

A charged track crossing a pixel sensor at an angle lights up several pixel columns in the local-x
(r-φ) direction **purely from geometry** — the more inclined the track, the more columns. A MIP
lights up roughly that many; a highly-ionizing MCP blooms charge into **extra** columns. So we
subtract the geometric expectation to isolate the bloom.

Per pixel hit (in `MCPProbQ::computeProbQ`):
1. Track direction in the module local frame: `lv = geomDet->toLocal(GlobalVector(px,py,pz))`.
   Local **z** is the sensor normal, local **x** is the r-φ direction (~100 µm pitch).
2. Geometric projected travel across the sensor in x: `Δx = thickness · (lv.x / lv.z)` = `thickness · tan(α_x)`.
3. **Predicted r-φ size** = `thickness · |lv.x / lv.z| / pitchX` (number of pixel columns from angle alone),
   with `thickness = geomDet->surface().bounds().thickness()` and
   `pitchX = PixelGeomDetUnit::specificTopology().pitch().first`.
4. **`pixSizeXresidual` = mean( measured `pixelCluster->sizeX()` − predicted )** over the track's pixel hits
   (`pixSizeXpred` stores the mean predicted size).

What "predicted" deliberately omits:
- **No Lorentz-drift / diffusion term.** The B field widens every cluster by a roughly constant
  ~1 pixel (Lorentz angle ~23° in BPIX); this is left out, so it appears as a constant baseline in
  the residual (QCD sits at ~0.8–1); since it shifts both signal and background equally, only the bloom discriminates.
- Track direction is the **global production momentum** (good near the beamline; mild approximation
  for low-momentum, i.e. high-Q reco, tracks).

`pixSizeXresidual` separates from QCD with AUC ≈ 0.87–0.88, **improving with charge** (the bloom keeps
growing while charge-based variables saturate), and is largely independent of `probQ_pixel`
(corr ≈ 0.15 in QCD) — a complementary, charge-robust second axis for a pass/fail design. It is not a
primary discriminant (dE/dx / probQ reach AUC ≈ 0.97).

### Other branches
- Track: `pt/eta/phi`, `ptError`, `normChi2`, `validFraction`, `charge`, `nValidPixelHits`,
  `nTrackerLayers`, `highPurity`, builtin `dedxPixel_builtin`/`dedxStrip_builtin`.
- Per-subdetector hit counts `nBPIX..nTEC`, `nPix`, `nStrip`; per-hit dE/dx vectors
  `pixelDedxHits`/`stripDedxHits`.
- Per-hit probXY diagnostics (`nPixClusters`, `nPixSpecInCPE`, `nPixXYvalid`, `pixXYrawMin`, …).
- Calorimeter (mostly not useful — see studies): `caloEmEnergy/HadEnergy` (matched calo-jet, an
  isolation proxy), packed-cand `pcCaloFrac/pcHcalFrac`.
- Gen match: `genMatched`, `gen_pt/eta/charge/pdgId`, `dR`, `chargeFromCurvature = gen_pt / reco_pt`
  (peaks at the true charge since reco assumes |q|=1). The `gen` tree holds one entry per gen MCP.

## Test / study scripts (`test/`)
These use PyROOT, so they need the CMSSW ROOT on the path — run `cmsenv` from
`CMSSW_15_0_19_patch2/src` first (i.e. `cd CMSSW_15_0_19_patch2/src && cmsenv`), then invoke the
script from `~/HSCP-MCP/`.
- `cmsstyle.py` — shared square CMS-style plotting helper (`apply`, `square_canvas`, `cms_label`, `save` → png+pdf).
- `plotCharge.py` — `chargeFromCurvature` (charge-scaling validation). Overlays one normalized
  `gen_pt/reco_pt` histogram per input file, each with a dashed line at its expected charge; the
  x-range auto-extends to cover the largest expected charge. Usage:
  `python3 plotCharge.py <out_noext> "file.root|label|Qexp" ...`, e.g. the M=2000 ladder
  `"mcp_M2000_Q6_v2.root|Q=2e|2" "...Q30...|Q=10e|10" "...Q45...|Q=15e|15" "...Q72...|Q=24e|24" "...Q90...|Q=30e|30"`
  → `mcp_chargeFromCurvature_v2`.
- `simSystemCoverage.py` — GEN-SIM detector reach: which subsystems (tracker, ECAL, HCAL, muon) the MCP
  reaches, swept across the full charge ladder; produces a heatmap and reach-fraction-vs-charge summary.
- `plotSigVsBkg.py` — signal-vs-background 1D overlays and the 2D `probQ_pixel`-vs-`sizeXresidual`
  plane. The per-charge 1D plots (`sigVsBkg_<var>_Q*`) cover `probQ_strip`, `ih_pixel`, `ih_full`.
  For `probQ_pixel` and `sizeXresidual` it instead writes a single multi-charge overlay across the
  full M=2000 ladder (`sigVsBkg_probQpixel_overlay`, `sigVsBkg_sizeXresidual_overlay`): one signal
  histogram per charge (colour = charge, Q=2/8/10/15/24/30e) over the QCD background (dashed black),
  all normalized. A matching overlay for the **isolated-tracks collection size per event**
  (`sigVsBkg_nIsoTracks_overlay[_logy]`) is produced the same way — events (normalized) vs number of
  isolated tracks in the event, counting all tracks in the ntuple regardless of dE/dx selection.
  Each 2D plane is emitted as a
  per-sample COLZ (`plane_<tag>_Q*`, `plane_<tag>_<bkg>`) and a per-track scatter overlay
  (`..._scatter_Q*`); the COLZ plots include the Pearson correlation factor `ρ` as an on-canvas label.
  It also writes a ROC for the plane (`roc_probQpixel_sizeXresidual_Q*`):
  signal-eff vs background-eff for each axis (higher = signal-like, keep `v>=thr`), with the AUC and
  the optimal single-axis cut (max Youden's J) marked, plus the best 2D rectangular cut from a grid
  scan. The in-plane signal acceptance (fraction with `probQ_pixel` defined) is printed and labelled,
  since the ROC efficiencies are conditional on entering the plane. It additionally auto-discovers the
  M=2000 charge ladder in the output dir (`mcp_M2000_Q{24,30,45,72,90}_v2.root` = Q=8/10/15/24/30e)
  and overlays them: `roc_charges_byvariable` (probQ_pixel solid + sizeXresidual dashed, colour=charge)
  and `roc_charges_combined` (best `probQ_pixel`x`sizeXresidual` rectangular-cut envelope per charge).
  Usage:
  ```
  cd $CMSSW_BASE/.. && cmsenv
  python3 MCPAnalyzer/Analyzer/test/plotSigVsBkg.py \
      $HOME/HSCP-MCP <Qe> $HOME/HSCP-MCP/mcp_M2000_Q<tag>_v2.root \
      "QCD 50-80=$HOME/HSCP-MCP/mcp_QCD_Pt50to80_v2.root"
  ```
  e.g. for Q=10e (file tag Q30):
  ```
  python3 MCPAnalyzer/Analyzer/test/plotSigVsBkg.py \
      $HOME/HSCP-MCP 10 $HOME/HSCP-MCP/mcp_M2000_Q30_v2.root \
      "QCD 50-80=$HOME/HSCP-MCP/mcp_QCD_Pt50to80_v2.root"
  ```
  The script auto-discovers all ladder files (`mcp_M2000_Q{6,24,30,45,72,90}_v2.root`) in the output
  dir for the multi-charge overlays and ROC comparisons; only the explicitly named signal file is used
  for the single-charge per-`Qe` plots.

## Status notes
- The calorimeter ionization deposit is large at SIM level (≲48 GeV for Q≲15e) but is **not
  accessible via standard MiniAOD objects** (it's energy along a charged track) — would need AOD.
