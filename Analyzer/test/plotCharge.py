#!/usr/bin/env python3
"""Charge-from-curvature plot (gen pT / reco pT) from MCPAnalyzer 'gen' trees.

Usage:
  python3 plotCharge.py <out_noext> "file1.root|label1|Qexp1" "file2.root|label2|Qexp2" ...
"""
import sys, ROOT
import cmsstyle

OUT = sys.argv[1]
specs = []
cols = [ROOT.kBlue + 1, ROOT.kRed + 1, ROOT.kGreen + 2, ROOT.kViolet + 1, ROOT.kOrange + 1, ROOT.kAzure + 1]
for i, a in enumerate(sys.argv[2:]):
    fn, lab, q = a.split("|")
    specs.append((lab, fn, int(q), cols[i % len(cols)]))

ROOT.gROOT.SetBatch(True)
cmsstyle.apply()
# x range covers the largest expected charge (peak sits at gen_pt/reco_pt ~ Q), ~0.1/bin
XHI = max(12.0, max(q for _, _, q, _ in specs) * 1.25)
NB = int(round(XHI / 0.1))
c = cmsstyle.square_canvas("c"); c.SetLeftMargin(0.13)
hs = ROOT.THStack("hs", ";charge from curvature  (gen p_{T} / reco p_{T});MCP tracks (normalized)")
leg = ROOT.TLegend(0.5, 0.72, 0.92, 0.88); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.030)
keep = []
lines = []
for lab, fn, qexp, col in specs:
    f = ROOT.TFile.Open(fn); gen = f.Get("MCPAnalyzer/gen")
    h = ROOT.TH1F("h_" + lab, "", NB, 0, XHI); h.SetDirectory(0)
    for e in gen:
        if e.matched and e.reco_pt > 0:
            h.Fill(e.gen_pt / e.reco_pt)
    if h.Integral() > 0: h.Scale(1. / h.Integral())
    h.SetLineColor(col); h.SetLineWidth(2)
    hs.Add(h); keep.append(h)
    leg.AddEntry(h, "%s  (expect %d)" % (lab, qexp), "l")
    f.Close()
hs.Draw("nostack hist")
hs.GetYaxis().SetTitleOffset(1.3)
ymax = hs.GetMaximum("nostack") * 1.15
for lab, fn, qexp, col in specs:
    l = ROOT.TLine(qexp, 0, qexp, ymax); l.SetLineColor(col); l.SetLineStyle(2); l.Draw(); lines.append(l)
leg.Draw()
cmsstyle.cms_label(c)
cmsstyle.save(c, OUT)
print("wrote", OUT + ".png/.pdf")
