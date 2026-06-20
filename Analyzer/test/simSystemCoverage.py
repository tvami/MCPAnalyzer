#!/usr/bin/env python3
"""Systematic SimHit-reach study vs charge at fixed mass.

For each charge point, find the MCP SimTracks (|type|==MCP_PDGID) and count their
PSimHits per subdetector. Builds:
  - a table of event-reach fraction (fraction of MCP events with >=1 MCP SimHit) per subdet
  - a CSV
  - a 2D heatmap (charge x subdetector) and a "deepest layer reached" curve vs charge

Usage:
  python3 systReach.py <MASS> <NFILES_PER_POINT> <outdir>
"""
import sys, glob, ROOT
from DataFormats.FWLite import Events, Handle

MASS    = int(sys.argv[1])
NFILES  = int(sys.argv[2])
OUTDIR  = sys.argv[3]
MCP_PDGID = 10000200
BASE = "/ceph/cms/store/user/tvami/HSCP_MCP"
CHARGES = [3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,60,72,90]

TK = [("BPIX","TrackerHitsPixelBarrelLowTof"),("BPIX","TrackerHitsPixelBarrelHighTof"),
      ("FPIX","TrackerHitsPixelEndcapLowTof"),("FPIX","TrackerHitsPixelEndcapHighTof"),
      ("TIB","TrackerHitsTIBLowTof"),("TIB","TrackerHitsTIBHighTof"),
      ("TID","TrackerHitsTIDLowTof"),("TID","TrackerHitsTIDHighTof"),
      ("TOB","TrackerHitsTOBLowTof"),("TOB","TrackerHitsTOBHighTof"),
      ("TEC","TrackerHitsTECLowTof"),("TEC","TrackerHitsTECHighTof")]
MU = [("DT","MuonDTHits"),("CSC","MuonCSCHits"),
      ("GEM","MuonGEMHits")]   # ME0 dropped: high-eta forward, acceptance-only (0 at all charges)
                               # RPC dropped: redundant with DT/CSC for the reach message
# Calorimeters use PCaloHit (geantTrackId, energy). "reached" = MCP deposits > CALO_ETHR GeV
# (an energy threshold, not >=1 hit, so the misleading front-face grazing at high Q is excluded).
CALO = [("ECAL","EcalHitsEB"),("ECAL","EcalHitsEE"),("HCAL","HcalHits")]
CALO_ETHR = 1.0
ALL = TK + MU
GROUPS = ["BPIX","FPIX","TIB","TID","TOB","TEC","ECAL","HCAL","DT","CSC","GEM"]

def find_files(q):
    pat = "%s/HSCP-MCP_M-%d-Q-%d_GEN-SIM_2024_v2_13p6TeV/*_GEN_*/*/0000/*GENSIM*.root" % (BASE, MASS, q)
    fs = sorted(glob.glob(pat))
    return ["file:" + f for f in fs[:NFILES]]

hsim = Handle("vector<SimTrack>")
hh   = {inst: Handle("vector<PSimHit>") for _, inst in ALL}
hc   = {inst: Handle("vector<PCaloHit>") for _, inst in CALO}

results = {}   # q -> dict(group->reach fraction)
hitavg  = {}   # q -> dict(group->avg MCP hits/event)
nevs    = {}
for q in CHARGES:
    fs = find_files(q)
    if not fs:
        print("Q-%d: no files" % q); continue
    tot = {g: 0 for g in GROUPS}; evc = {g: 0 for g in GROUPS}; nev = 0
    ev = Events(fs)
    for e in ev:
        e.getByLabel("g4SimHits", hsim)
        mcpIds = set(st.trackId() for st in hsim.product() if abs(st.type()) == MCP_PDGID)
        if not mcpIds:
            continue
        nev += 1; seen = set()
        for grp, inst in ALL:
            e.getByLabel("g4SimHits", inst, hh[inst])
            if not hh[inst].isValid():
                continue
            c = sum(1 for ph in hh[inst].product() if ph.trackId() in mcpIds)
            if c:
                tot[grp] += c; seen.add(grp)
        # calorimeters: sum MCP-deposited energy, "reached" if > threshold
        caloE = {"ECAL": 0.0, "HCAL": 0.0}
        for grp, inst in CALO:
            e.getByLabel("g4SimHits", inst, hc[inst])
            if not hc[inst].isValid():
                continue
            for ph in hc[inst].product():
                if ph.geantTrackId() in mcpIds:
                    caloE[grp] += ph.energy()
        for grp in ("ECAL", "HCAL"):
            if caloE[grp] > CALO_ETHR:
                seen.add(grp)
        for g in seen:
            evc[g] += 1
    results[q] = {g: (evc[g] / nev if nev else 0.) for g in GROUPS}
    hitavg[q]  = {g: (tot[g] / nev if nev else 0.) for g in GROUPS}
    nevs[q] = nev
    print("Q-%-3d (%2de) nev=%4d | " % (q, q // 3, nev) +
          " ".join("%s:%4.0f%%" % (g, 100. * results[q][g]) for g in GROUPS))

# CSV
csv = OUTDIR + "/systReach_M%d.csv" % MASS
with open(csv, "w") as fo:
    fo.write("charge_e,nev," + ",".join(GROUPS) + "\n")
    for q in CHARGES:
        if q not in results: continue
        fo.write("%d,%d," % (q // 3, nevs[q]) + ",".join("%.4f" % results[q][g] for g in GROUPS) + "\n")
print("wrote", csv)

# ---------- plots (square, CMS style, png+pdf) ----------
ROOT.gROOT.SetBatch(True)
import cmsstyle
cmsstyle.apply()
qs = [q for q in CHARGES if q in results]

# 2D heatmap: y=charge, x=subdet, z=event-reach fraction
ROOT.gStyle.SetPaintTextFormat(".0f")
h2 = ROOT.TH2F("reach", "", len(GROUPS), 0, len(GROUPS), len(qs), 0, len(qs))
for i, g in enumerate(GROUPS):
    h2.GetXaxis().SetBinLabel(i + 1, g)
for j, q in enumerate(qs):
    h2.GetYaxis().SetBinLabel(j + 1, "%d" % (q // 3))
    for i, g in enumerate(GROUPS):
        h2.SetBinContent(i + 1, j + 1, 100. * results[q][g])
h2.GetYaxis().SetTitle("charge Q [e]"); h2.GetZaxis().SetTitle("MCP events reaching subdet [%]")
h2.GetXaxis().SetLabelSize(0.038); h2.GetYaxis().SetLabelSize(0.034)
h2.GetZaxis().SetTitleOffset(1.2); h2.SetMarkerSize(1.0); h2.SetMinimum(0); h2.SetMaximum(100)
c = cmsstyle.square_canvas("c", right=0.15); c.SetLeftMargin(0.12)
h2.Draw("COLZ TEXT")
ln = ROOT.TLine(6, 0, 6, len(qs)); ln.SetLineColor(ROOT.kWhite); ln.SetLineWidth(2); ln.Draw()
cmsstyle.cms_label(c)
cmsstyle.save(c, OUTDIR + "/systReach_M%d_heatmap" % MASS)
print("wrote heatmap")

# summary line plot: key boundaries vs charge (legend in upper-right gap)
def colv(name): return [100. * results[q][name] for q in qs]
xq = [q // 3 for q in qs]
series = [("BPIX (pixel)", colv("BPIX"), ROOT.kBlack, 20),
          ("TOB (outer strip)", colv("TOB"), ROOT.kGreen + 2, 21),
          ("TEC (fwd strip)", colv("TEC"), ROOT.kOrange + 1, 22),
          ("ECAL (>1 GeV)", colv("ECAL"), ROOT.kCyan + 2, 47),
          ("HCAL (>1 GeV)", colv("HCAL"), ROOT.kMagenta + 2, 43),
          ("DT (muon barrel)", colv("DT"), ROOT.kRed + 1, 23)]
c2 = cmsstyle.square_canvas("c2"); c2.SetLeftMargin(0.13); c2.SetGridx(); c2.SetGridy()
fr = c2.DrawFrame(0, 0, max(xq) + 2, 112)
fr.GetXaxis().SetTitle("charge Q [e]"); fr.GetYaxis().SetTitle("fraction of MCP events reaching subdet [%]")
fr.GetYaxis().SetTitleOffset(1.3)
leg = ROOT.TLegend(0.61, 0.56, 0.93, 0.80); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.030)
graphs = []
for name, ys, cc, mk in series:
    g = ROOT.TGraph()
    for i, (q, y) in enumerate(zip(xq, ys)): g.SetPoint(i, q, y)
    g.SetLineColor(cc); g.SetMarkerColor(cc); g.SetMarkerStyle(mk); g.SetLineWidth(2); g.SetMarkerSize(1.2)
    g.Draw("PL same"); graphs.append(g); leg.AddEntry(g, name, "lp")
leg.Draw()
l50 = ROOT.TLine(0, 50, max(xq) + 2, 50); l50.SetLineStyle(2); l50.SetLineColor(ROOT.kGray + 2); l50.Draw()
cmsstyle.cms_label(c2)
cmsstyle.save(c2, OUTDIR + "/systReach_M%d_summary" % MASS)
print("wrote summary")
