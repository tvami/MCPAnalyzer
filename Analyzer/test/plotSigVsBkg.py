#!/usr/bin/env python3
"""Overlay signal (gen-matched MCP) vs background (all measured tracks) on the
discriminant variables, plus the 2D probQ_pixel-vs-sizeXresidual planes.

Usage:
  python3 plotSigVsBkg.py <outdir> <sigQe> <sig.root> "label1=bkg1.root" [...]
e.g.
  python3 plotSigVsBkg.py $HOME/HSCP-MCP 10 mcp_M2000_Q30_v2.root "QCD 50-80=mcp_QCD_Pt50to80_v2.root"

probQ_pixel / probQ_strip convention (after the "never combine a per-hit 0" fix):
  p <  0  -> UNDEFINED (no usable charge hit)  -> excluded from every plot
  p == 0  -> genuine float underflow           -> placed in the last visible -log10 bin
  else    -> -log10(p), clamped to CAP
"""
import sys, os, math, ROOT
 # find cmsstyle when run elsewhere
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))) 
import cmsstyle
ROOT.gROOT.SetBatch(True); cmsstyle.apply()

OUT = sys.argv[1]
SIGQ = int(sys.argv[2])
SIGF = sys.argv[3]
BKGS = []
for a in sys.argv[4:]:
    lab, fn = a.split("=", 1); BKGS.append((lab, fn))

CAP = 15.0
def nlog(p):
    """-log10(probQ) for plotting; None if undefined (p<0), CAP for genuine underflow (p==0)."""
    if p < 0:   return None
    if p == 0:  return CAP - 0.05
    return max(0.0, min(CAP - 0.05, -math.log10(p)))

# selection: measured pixel tracks. signal also requires gen-match.
def sig_ok(e):  return e.genMatched and e.hasDeDx and e.nPix > 0
def bkg_ok(e):  return e.hasDeDx and e.nPix > 0

SIGCOL = ROOT.kRed + 1
BKGCOLS = [ROOT.kBlack, ROOT.kAzure + 1, ROOT.kGreen + 2, ROOT.kOrange + 1]

def load_event_niso(fn):
    """Distribution of isolated-track multiplicity per event (all TTree entries, no selection)."""
    f = ROOT.TFile.Open(fn); t = f.Get("MCPAnalyzer/tracks")
    counts = {}
    for e in t:
        key = (int(e.run), int(e.lumi), int(e.event))
        counts[key] = counts.get(key, 0) + 1
    h = ROOT.TH1F("h_niso", "", 50, 0.5, 50.5); h.SetDirectory(0)
    for v in counts.values():
        h.Fill(v)
    if h.Integral() > 0: h.Scale(1. / h.Integral())
    f.Close(); return h

def load1d(fn, sel, getter, nb, lo, hi):
    f = ROOT.TFile.Open(fn); t = f.Get("MCPAnalyzer/tracks")
    h = ROOT.TH1F("h", "", nb, lo, hi); h.SetDirectory(0)
    for e in t:
        if not sel(e): continue
        v = getter(e)
        if v is None: continue
        h.Fill(v)
    if h.Integral() > 0: h.Scale(1. / h.Integral())
    f.Close(); return h

# ---- 2D planes (COLZ per-sample + contour overlay + scatter overlay), one styling for all ----
# Each plane: (tag, xtitle, xget, (ylo,yhi), ytitle, yget). x/y getters return None to skip a track.
PLANES = [
    ("probQpixel_sizeXresidual",
     "-log_{10}(probQ_{pixel})", lambda e: nlog(e.probQ_pixel),
     (-2, 8), "pixel sizeX residual [pixels]",
     lambda e: (None if e.pixSizeXresidual < -90 else e.pixSizeXresidual)),
]

def load2d(fn, sel, xget, yget, ylo, yhi):
    f = ROOT.TFile.Open(fn); t = f.Get("MCPAnalyzer/tracks")
    h = ROOT.TH2F("h2", "", 40, 0, CAP, 40, ylo, yhi); h.SetDirectory(0)
    g = ROOT.TGraph()
    for e in t:
        if not sel(e): continue
        x, y = xget(e), yget(e)
        if x is None or y is None: continue
        h.Fill(x, y); g.SetPoint(g.GetN(), x, y)
    f.Close(); return h, g

def draw_colz(h, xt, yt, title, out, corr=None):
    c = cmsstyle.square_canvas("c", right=0.15); c.SetLeftMargin(0.13)
    h.GetXaxis().SetTitle(xt); h.GetYaxis().SetTitle(yt)
    h.GetZaxis().SetTitle("Tracks"); h.GetYaxis().SetTitleOffset(1.3)
    h.Draw("COLZ"); cmsstyle.cms_label(c)
    lt = ROOT.TLatex(); lt.SetNDC(); lt.SetTextSize(0.035); lt.DrawLatex(0.16, 0.86, title)
    if corr is not None:
        lt.SetTextSize(0.030); lt.DrawLatex(0.16, 0.80, "#rho = %.3f" % corr)
    cmsstyle.save(c, out); print("wrote", os.path.basename(out), "entries=%d" % int(h.GetEntries()))

for tag, xt, xget, (ylo, yhi), yt, yget in PLANES:
    hsig, gsig = load2d(SIGF, sig_ok, xget, yget, ylo, yhi)
    bkg2d = [(lab, load2d(fn, bkg_ok, xget, yget, ylo, yhi)) for lab, fn in BKGS]

    # per-sample COLZ
    draw_colz(hsig, xt, yt, "M=2000, Q=%de signal" % SIGQ, "%s/plane_%s_Q%d" % (OUT, tag, SIGQ),
              corr=hsig.GetCorrelationFactor())
    for lab, (h, _) in bkg2d:
        draw_colz(h, xt, yt, lab, "%s/plane_%s_%s" % (OUT, tag, lab.split()[0]),
                  corr=h.GetCorrelationFactor())

    # scatter overlay (background first, signal on top)
    c = cmsstyle.square_canvas("c"); c.SetLeftMargin(0.13)
    fr = c.DrawFrame(0, ylo, CAP, yhi)
    fr.GetXaxis().SetTitle(xt); fr.GetYaxis().SetTitle(yt); fr.GetYaxis().SetTitleOffset(1.3)
    gk = []
    for i, (lab, (_, g)) in enumerate(bkg2d):
        g.SetMarkerColor(BKGCOLS[i % len(BKGCOLS)]); g.SetMarkerStyle(24); g.SetMarkerSize(0.6)
        g.Draw("P same"); gk.append((lab, g))
    gsig.SetMarkerColor(ROOT.kRed + 1); gsig.SetMarkerStyle(20); gsig.SetMarkerSize(0.6); gsig.Draw("P same")
    leg = ROOT.TLegend(0.5, 0.72, 0.92, 0.88); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.028)
    leg.AddEntry(gsig, "signal Q=%de (MCP)  [%d]" % (SIGQ, gsig.GetN()), "p")
    for lab, g in gk: leg.AddEntry(g, "%s  [%d]" % (lab, g.GetN()), "p")
    leg.Draw(); cmsstyle.cms_label(c)
    cmsstyle.save(c, "%s/plane_%s_scatter_Q%d" % (OUT, tag, SIGQ))
    print("wrote plane_%s_scatter_Q%d" % (tag, SIGQ))

# ============================================================================
# ROC for the probQ_pixel-vs-sizeXresidual plane: best cut on each axis + best 2D
# rectangular cut. Both variables are "higher = more signal-like", so a cut keeps
# v >= thr. Background = first sample given on the command line.
# Population = tracks where BOTH variables are defined (same as the scatter plane);
# the fraction of signal that even enters this plane (acceptance) is reported too.
# ============================================================================
import bisect
BKGLAB, BKGFN = BKGS[0]

def collect_pairs(fn, sel):
    f = ROOT.TFile.Open(fn); t = f.Get("MCPAnalyzer/tracks")
    pairs, ntot = [], 0
    for e in t:
        if not sel(e): continue
        ntot += 1
        x = nlog(e.probQ_pixel)
        y = None if e.pixSizeXresidual < -90 else e.pixSizeXresidual
        if x is None or y is None: continue
        pairs.append((x, y))
    f.Close(); return pairs, ntot

sigp, sig_ntot = collect_pairs(SIGF, sig_ok)
bkgp, bkg_ntot = collect_pairs(BKGFN, bkg_ok)
print("ROC population: signal %d/%d in-plane (acc=%.2f), bkg %d/%d (acc=%.2f)"
      % (len(sigp), sig_ntot, len(sigp) / max(sig_ntot, 1),
         len(bkgp), bkg_ntot, len(bkgp) / max(bkg_ntot, 1)))

def roc1d(svals, bvals):
    ns, nb = len(svals), len(bvals)
    ss, bs = sorted(svals), sorted(bvals)
    pts = [(1.0, 1.0, float("-inf"))]
    for thr in sorted(set(svals + bvals)):
        se = (ns - bisect.bisect_left(ss, thr)) / ns
        be = (nb - bisect.bisect_left(bs, thr)) / nb
        pts.append((be, se, thr))
    pts.append((0.0, 0.0, float("inf")))
    return pts

def auc_of(pts):
    p = sorted((be, se) for be, se, _ in pts)
    a = 0.0
    for i in range(1, len(p)):
        a += (p[i][0] - p[i - 1][0]) * (p[i][1] + p[i - 1][1]) / 2.0
    return a

def best_J(pts):  # max Youden's J = sigEff - bkgEff
    return max(pts, key=lambda t: t[1] - t[0])  # (be, se, thr)

AXES = [("-log_{10}(probQ_{pixel})", "probQpixel", 0, [p[0] for p in sigp], [p[0] for p in bkgp], ROOT.kRed + 1, 20),
        ("pixel sizeX residual",       "sizeXresidual", 1, [p[1] for p in sigp], [p[1] for p in bkgp], ROOT.kAzure + 1, 21)]

c = cmsstyle.square_canvas("c"); c.SetLeftMargin(0.13); c.SetGridx(); c.SetGridy()
fr = c.DrawFrame(0, 0, 1, 1.02)
fr.GetXaxis().SetTitle("Background (%s) efficiency" % BKGLAB)
fr.GetYaxis().SetTitle("Signal efficiency"); fr.GetYaxis().SetTitleOffset(1.3)
diag = ROOT.TLine(0, 0, 1, 1); diag.SetLineStyle(2); diag.SetLineColor(ROOT.kGray + 2); diag.Draw()
leg = ROOT.TLegend(0.4, 0.18, 0.92, 0.40); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.027)
keep = []
best = {}
for title, key, idx, sv, bv, col, mk in AXES:
    pts = roc1d(sv, bv)
    a = auc_of(pts)
    g = ROOT.TGraph()
    for be, se, _ in sorted((p[0], p[1], p[2]) for p in pts):
        g.SetPoint(g.GetN(), be, se)
    g.SetLineColor(col); g.SetLineWidth(3); g.SetMarkerColor(col); g.Draw("L same"); keep.append(g)
    be, se, thr = best_J(pts); best[key] = (be, se, thr)
    m = ROOT.TMarker(be, se, mk); m.SetMarkerColor(col); m.SetMarkerSize(1.6); m.Draw(); keep.append(m)
    leg.AddEntry(g, "%s  (AUC=%.3f)" % (title, a), "l")
    print("  %-22s AUC=%.3f  best cut: keep v>=%.2f  -> sigEff=%.3f bkgEff=%.3f (rej=%.3f)"
          % (title, a, thr, se, be, 1 - be))
leg.Draw(); cmsstyle.cms_label(c)

# 2D rectangular cut scan: keep x>=xt AND y>=yt, maximize Youden's J
xg = [i * 0.5 for i in range(0, int(2 * CAP) + 1)]
yg = [i * 0.25 for i in range(-4, 25)]
ns, nb = len(sigp), len(bkgp)
b2 = (-1, None, None, None, None)
for xt in xg:
    for yt in yg:
        se = sum(1 for x, y in sigp if x >= xt and y >= yt) / ns
        be = sum(1 for x, y in bkgp if x >= xt and y >= yt) / nb
        J = se - be
        if J > b2[0]: b2 = (J, xt, yt, se, be)
_, bxt, byt, bse, bbe = b2
print("  2D best cut: -log10(probQ_pixel)>=%.2f AND sizeXresidual>=%.2f -> sigEff=%.3f bkgEff=%.3f (rej=%.3f)"
      % (bxt, byt, bse, bbe, 1 - bbe))
m2 = ROOT.TMarker(bbe, bse, 29); m2.SetMarkerColor(ROOT.kGreen + 2); m2.SetMarkerSize(2.4); m2.Draw(); keep.append(m2)
leg.AddEntry(m2, "2D cut: probQpix>=%.1f & szX>=%.2f" % (bxt, byt), "p")
lt = ROOT.TLatex(); lt.SetNDC(); lt.SetTextSize(0.026)
lt.DrawLatex(0.16, 0.86, "signal Q=%de   (in-plane acc %.0f%%)" % (SIGQ, 100 * len(sigp) / max(sig_ntot, 1)))
cmsstyle.save(c, "%s/roc_probQpixel_sizeXresidual_Q%d" % (OUT, SIGQ))
print("wrote roc_probQpixel_sizeXresidual_Q%d" % SIGQ)

# ============================================================================
# Multi-charge ROC overlay (same plane, same QCD background). Auto-discovers the
# standard M=2000 ntuples in OUT.  Q[e] -> file tag is 3*Q (8e->Q24 ... 30e->Q90).
#   (1) per-variable overlay: probQ_pixel SOLID, sizeXresidual DASHED, colour = charge
#   (2) combined overlay: optimal axis-aligned 2D-cut envelope, one curve per charge
# ============================================================================
CHARGES = [(8, 24), (10, 30), (15, 45), (24, 72), (30, 90)]
QCOL = {8: ROOT.kViolet + 1, 10: ROOT.kAzure + 1, 15: ROOT.kGreen + 2,
        24: ROOT.kOrange + 1, 30: ROOT.kRed + 1}

charge_pairs = []  # (Qe, sigpairs, acc)
for Qe, tag in CHARGES:
    fn = os.path.join(OUT, "mcp_M2000_Q%d_v2.root" % tag)
    if not os.path.exists(fn):
        print("  (skip Q=%de: %s missing)" % (Qe, os.path.basename(fn))); continue
    sp, nt = collect_pairs(fn, sig_ok)
    if sp: charge_pairs.append((Qe, sp, len(sp) / max(nt, 1)))

def trap_auc(pts):  # pts = list of (be, se) sorted by be
    a = 0.0
    for i in range(1, len(pts)):
        a += (pts[i][0] - pts[i - 1][0]) * (pts[i][1] + pts[i - 1][1]) / 2.0
    return a

# ---- (1) per-variable overlay ----
c = cmsstyle.square_canvas("c"); c.SetLeftMargin(0.13); c.SetGridx(); c.SetGridy()
fr = c.DrawFrame(0, 0, 1, 1.02)
fr.GetXaxis().SetTitle("Background (%s) efficiency" % BKGLAB)
fr.GetYaxis().SetTitle("Signal efficiency"); fr.GetYaxis().SetTitleOffset(1.3)
ROOT.TLine(0, 0, 1, 1).DrawClone("").SetLineStyle(2)
keep = []
legQ = ROOT.TLegend(0.55, 0.13, 0.92, 0.13 + 0.05 * len(charge_pairs))
legQ.SetBorderSize(0); legQ.SetFillStyle(0); legQ.SetTextSize(0.027)
for Qe, sp, acc in charge_pairs:
    bvx = [p[0] for p in bkgp]; bvy = [p[1] for p in bkgp]
    for idx, style in ((0, 1), (1, 2)):  # x=probQ solid, y=sizeX dashed
        pts = roc1d([p[idx] for p in sp], (bvx if idx == 0 else bvy))
        g = ROOT.TGraph()
        for be, se, _ in sorted((q[0], q[1], q[2]) for q in pts): g.SetPoint(g.GetN(), be, se)
        g.SetLineColor(QCOL[Qe]); g.SetLineStyle(style); g.SetLineWidth(3)
        g.Draw("L same"); keep.append(g)
        if idx == 0:
            legQ.AddEntry(g, "Q=%de (acc %.0f%%)" % (Qe, 100 * acc), "l")
            print("  Q=%-2de probQ_pixel AUC=%.3f  sizeXresidual AUC=%.3f"
                  % (Qe, trap_auc(sorted((p[0], p[1]) for p in roc1d([q[0] for q in sp], bvx))),
                         trap_auc(sorted((p[0], p[1]) for p in roc1d([q[1] for q in sp], bvy)))))
legQ.Draw()
# style key (solid vs dashed)
ls = ROOT.TLegend(0.55, 0.40, 0.92, 0.50); ls.SetBorderSize(0); ls.SetFillStyle(0); ls.SetTextSize(0.026)
g1 = ROOT.TGraph(); g1.SetLineStyle(1); g1.SetLineWidth(3); g1.SetLineColor(ROOT.kBlack)
g2 = ROOT.TGraph(); g2.SetLineStyle(2); g2.SetLineWidth(3); g2.SetLineColor(ROOT.kBlack)
keep += [g1, g2]; ls.AddEntry(g1, "probQ_{pixel} (solid)", "l"); ls.AddEntry(g2, "sizeX residual (dashed)", "l"); ls.Draw()
cmsstyle.cms_label(c); cmsstyle.save(c, "%s/roc_charges_byvariable" % OUT)
print("wrote roc_charges_byvariable")

# ---- (2) combined (2D-cut envelope) overlay, log-x, with the optimum cut marked ----
def grid_frac(pairs):
    n = len(pairs) or 1
    return {(xt, yt): sum(1 for x, y in pairs if x >= xt and y >= yt) / n for xt in xg for yt in yg}
bkg_grid = grid_frac(bkgp)
XMIN = 1e-3   # log-x floor (smallest resolvable bkg eff ~ 1/Nbkg)

c = cmsstyle.square_canvas("c"); c.SetLeftMargin(0.13); c.SetGridx(); c.SetGridy(); c.SetLogx()
fr = c.DrawFrame(XMIN, 0, 1, 1.02)
fr.GetXaxis().SetTitle("Background (%s) efficiency" % BKGLAB)
fr.GetYaxis().SetTitle("Signal efficiency"); fr.GetYaxis().SetTitleOffset(1.3)
# X axis: decade labels only (dropping SetMoreLogLabels removes the crowded 2x10^n labels)
keep2 = []
legC = ROOT.TLegend(0.15, 0.14, 0.80, 0.14 + 0.052 * len(charge_pairs))
legC.SetBorderSize(0); legC.SetFillStyle(0); legC.SetTextSize(0.023)
for Qe, sp, acc in charge_pairs:
    sg = grid_frac(sp)
    raw = sorted(set((bkg_grid[k], sg[k]) for k in sg) | {(XMIN, 0.), (1., 1.)})
    env, bestse = [], -1.0
    for be, se in raw:
        bestse = max(bestse, se); env.append((max(be, XMIN), bestse))
    g = ROOT.TGraph()
    for be, se in env: g.SetPoint(g.GetN(), be, se)
    g.SetLineColor(QCOL[Qe]); g.SetLineWidth(3); g.Draw("L same"); keep2.append(g)
    # optimum = grid cut maximizing Youden's J = sigEff - bkgEff
    bk = max(sg, key=lambda k: sg[k] - bkg_grid[k])
    bse, bbe = sg[bk], bkg_grid[bk]
    star = ROOT.TMarker(max(bbe, XMIN), bse, 29); star.SetMarkerColor(QCOL[Qe]); star.SetMarkerSize(2.6)
    star.Draw(); keep2.append(star)
    legC.AddEntry(g, "Q=%de: pQ#geq%.1f, szX#geq%.2f (#varepsilon_{s}=%.2f, #varepsilon_{b}=%.3f, AUC=%.3f)"
                  % (Qe, bk[0], bk[1], bse, bbe, trap_auc(env)), "l")
    print("  Q=%-2de combined opt: probQpix>=%.2f sizeX>=%.2f -> sigEff=%.3f bkgEff=%.3f (rej=%.3f) AUC=%.3f"
          % (Qe, bk[0], bk[1], bse, bbe, 1 - bbe, trap_auc(env)))
legC.Draw()
lt = ROOT.TLatex(); lt.SetNDC(); lt.SetTextSize(0.024)
lt.DrawLatex(0.15, 0.46, "combined = best probQ_{pixel}#timessizeX cut;  #star = optimum")
cmsstyle.cms_label(c); cmsstyle.save(c, "%s/roc_charges_combined" % OUT)
print("wrote roc_charges_combined")

# ---- (3) sizeXresidual-only overlay, log-x, with the optimum cut marked ----
bvy = [p[1] for p in bkgp]
c = cmsstyle.square_canvas("c"); c.SetLeftMargin(0.13); c.SetGridx(); c.SetGridy(); c.SetLogx()
fr = c.DrawFrame(XMIN, 0, 1, 1.02)
fr.GetXaxis().SetTitle("Background (%s) efficiency" % BKGLAB)
fr.GetYaxis().SetTitle("Signal efficiency"); fr.GetYaxis().SetTitleOffset(1.3)
keep3 = []
legS = ROOT.TLegend(0.15, 0.14, 0.80, 0.14 + 0.052 * len(charge_pairs))
legS.SetBorderSize(0); legS.SetFillStyle(0); legS.SetTextSize(0.024)
for Qe, sp, acc in charge_pairs:
    pts = roc1d([p[1] for p in sp], bvy)
    g = ROOT.TGraph()
    for be, se, _ in sorted((q[0], q[1], q[2]) for q in pts): g.SetPoint(g.GetN(), max(be, XMIN), se)
    g.SetLineColor(QCOL[Qe]); g.SetLineWidth(3); g.Draw("L same"); keep3.append(g)
    bbe, bse, bthr = best_J(pts)
    star = ROOT.TMarker(max(bbe, XMIN), bse, 29); star.SetMarkerColor(QCOL[Qe]); star.SetMarkerSize(2.6)
    star.Draw(); keep3.append(star)
    legS.AddEntry(g, "Q=%de: szX#geq%.2f (#varepsilon_{s}=%.2f, #varepsilon_{b}=%.3f, AUC=%.3f)"
                  % (Qe, bthr, bse, bbe, trap_auc(sorted((q[0], q[1]) for q in pts))), "l")
    print("  Q=%-2de sizeX opt: sizeX>=%.2f -> sigEff=%.3f bkgEff=%.3f (rej=%.3f) AUC=%.3f"
          % (Qe, bthr, bse, bbe, 1 - bbe, trap_auc(sorted((q[0], q[1]) for q in pts))))
legS.Draw()
lt = ROOT.TLatex(); lt.SetNDC(); lt.SetTextSize(0.024)
lt.DrawLatex(0.15, 0.46, "pixel sizeX residual only;  #star = optimum")
cmsstyle.cms_label(c); cmsstyle.save(c, "%s/roc_charges_sizeX" % OUT)
print("wrote roc_charges_sizeX")

# ============================================================================
# Multi-charge SIGNAL overlays: one plot per variable showing the full M=2000
# charge ladder (colour = charge) over the QCD background (dashed black), all
# normalized. Replaces the per-charge sigVsBkg_<var>_Q* for these variables.
# Auto-discovers mcp_M2000_Q{tag}_v2.root in OUT (tag = 3 x charge[e]).
# ============================================================================
OVL = [(2, 6), (8, 24), (10, 30), (15, 45), (24, 72), (30, 90)]
OVLCOL = {2: ROOT.kBlue + 1, 8: ROOT.kViolet + 1, 10: ROOT.kAzure + 1,
          15: ROOT.kGreen + 2, 24: ROOT.kOrange + 1, 30: ROOT.kRed + 1}
OVL_VARS = [
    ("probQpixel", 45, 0, CAP, "-log_{10}(probQ_{pixel})   (#rightarrow signal-like)",
     lambda e: nlog(e.probQ_pixel)),
    ("sizeXresidual", 40, -2, 8, "pixel sizeX residual [pixels]",
     lambda e: (None if e.pixSizeXresidual < -90 else e.pixSizeXresidual)),
]
for key, nb, lo, hi, xt, get in OVL_VARS:
    c = cmsstyle.square_canvas("c"); c.SetLeftMargin(0.13)
    hbkg = load1d(BKGFN, bkg_ok, get, nb, lo, hi)
    hbkg.SetLineColor(ROOT.kBlack); hbkg.SetLineWidth(2); hbkg.SetLineStyle(2)
    hsigs = []
    for Qe, tag in OVL:
        fn = os.path.join(OUT, "mcp_M2000_Q%d_v2.root" % tag)
        if not os.path.exists(fn):
            print("  (overlay skip Q=%de: %s missing)" % (Qe, os.path.basename(fn))); continue
        h = load1d(fn, sig_ok, get, nb, lo, hi)
        h.SetLineColor(OVLCOL[Qe]); h.SetLineWidth(3)
        hsigs.append((Qe, h))
    hmax = max([hbkg.GetMaximum()] + [h.GetMaximum() for _, h in hsigs])
    # smallest non-zero bin across all curves -> sensible logY floor
    hmin = min([h.GetBinContent(b) for h in [hbkg] + [hh for _, hh in hsigs]
                for b in range(1, h.GetNbinsX() + 1) if h.GetBinContent(b) > 0] or [1e-4])
    frame = hsigs[0][1] if hsigs else hbkg
    for logy in (False, True):
        c = cmsstyle.square_canvas("c"); c.SetLeftMargin(0.13); c.SetLogy(logy)
        frame.SetMaximum(hmax * (5.0 if logy else 1.25))
        frame.SetMinimum(hmin / 2.0 if logy else 0)
        frame.GetXaxis().SetTitle(xt); frame.GetYaxis().SetTitle("Tracks (normalized)")
        frame.GetYaxis().SetTitleOffset(1.3)
        frame.Draw("hist")
        for _, h in hsigs: h.Draw("hist same")
        hbkg.Draw("hist same")
        leg = ROOT.TLegend(0.5, 0.60, 0.92, 0.88); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.026)
        for Qe, h in hsigs: leg.AddEntry(h, "signal Q=%de (MCP)" % Qe, "l")
        leg.AddEntry(hbkg, BKGLAB, "l")
        leg.Draw(); cmsstyle.cms_label(c)
        out = "sigVsBkg_%s_overlay%s" % (key, "_logy" if logy else "")
        cmsstyle.save(c, "%s/%s" % (OUT, out)); print("wrote", out)

# ---- Isolated tracks collection size overlay ----
hbkg_n = load_event_niso(BKGFN)
hbkg_n.SetLineColor(ROOT.kBlack); hbkg_n.SetLineWidth(2); hbkg_n.SetLineStyle(2)
hsigs_n = []
for Qe, tag in OVL:
    fn = os.path.join(OUT, "mcp_M2000_Q%d_v2.root" % tag)
    if not os.path.exists(fn):
        print("  (niso overlay skip Q=%de: %s missing)" % (Qe, os.path.basename(fn))); continue
    h = load_event_niso(fn)
    h.SetLineColor(OVLCOL[Qe]); h.SetLineWidth(3)
    hsigs_n.append((Qe, h))
if hsigs_n or True:
    hmax_n = max([hbkg_n.GetMaximum()] + [h.GetMaximum() for _, h in hsigs_n])
    hmin_n = min([h.GetBinContent(b) for h in [hbkg_n] + [hh for _, hh in hsigs_n]
                  for b in range(1, h.GetNbinsX() + 1) if h.GetBinContent(b) > 0] or [1e-4])
    frame_n = hsigs_n[0][1] if hsigs_n else hbkg_n
    for logy in (False, True):
        c = cmsstyle.square_canvas("c"); c.SetLeftMargin(0.13); c.SetLogy(logy)
        frame_n.SetMaximum(hmax_n * (5.0 if logy else 1.25))
        frame_n.SetMinimum(hmin_n / 2.0 if logy else 0)
        frame_n.GetXaxis().SetTitle("Isolated tracks per event")
        frame_n.GetYaxis().SetTitle("Events (normalized)")
        frame_n.GetYaxis().SetTitleOffset(1.3)
        frame_n.Draw("hist")
        for _, h in hsigs_n: h.Draw("hist same")
        hbkg_n.Draw("hist same")
        leg = ROOT.TLegend(0.5, 0.60, 0.92, 0.88); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.026)
        for Qe, h in hsigs_n: leg.AddEntry(h, "signal Q=%de (MCP)" % Qe, "l")
        leg.AddEntry(hbkg_n, BKGLAB, "l")
        leg.Draw(); cmsstyle.cms_label(c)
        out = "sigVsBkg_nIsoTracks_overlay%s" % ("_logy" if logy else "")
        cmsstyle.save(c, "%s/%s" % (OUT, out)); print("wrote", out)
