"""Minimal shared CMS plotting style for the MCPAnalyzer test plots.

Usage:
    import cmsstyle
    cmsstyle.apply()
    c = cmsstyle.square_canvas("c")
    ...
    cmsstyle.cms_label(c)
    cmsstyle.save(c, "/path/plotname")   # writes .png and .pdf
"""
import ROOT


def apply():
    s = ROOT.gStyle
    s.SetOptStat(0)
    s.SetOptTitle(0)
    s.SetCanvasColor(0)
    s.SetPadColor(0)
    s.SetFrameBorderMode(0)
    s.SetCanvasBorderMode(0)
    s.SetPadTickX(1)
    s.SetPadTickY(1)
    s.SetLabelFont(42, "xyz")
    s.SetTitleFont(42, "xyz")
    s.SetTextFont(42)
    s.SetLabelSize(0.04, "xyz")
    s.SetTitleSize(0.045, "xyz")
    s.SetNdivisions(510, "xyz")
    s.SetPalette(ROOT.kBird)


_canvas_count = [0]
def square_canvas(name="c", size=820, right=0.05):
    # auto-uniquify the ROOT name so repeated square_canvas("c") calls don't
    # trigger "Deleting canvas with same name" warnings
    _canvas_count[0] += 1
    c = ROOT.TCanvas("%s_%d" % (name, _canvas_count[0]), "", size, size)
    c.SetLeftMargin(0.14)
    c.SetRightMargin(right)
    c.SetTopMargin(0.08)
    c.SetBottomMargin(0.13)
    return c


def cms_label(pad, sim=True, energy="13.6 TeV"):
    pad.cd()
    top = 1 - pad.GetTopMargin()
    x = pad.GetLeftMargin()
    l = ROOT.TLatex(); l.SetNDC(); l.SetTextFont(61); l.SetTextSize(0.05)
    l.DrawLatex(x, top + 0.012, "CMS")
    if sim:
        l2 = ROOT.TLatex(); l2.SetNDC(); l2.SetTextFont(52); l2.SetTextSize(0.038)
        l2.DrawLatex(x + 0.10, top + 0.012, "Simulation")
    r = ROOT.TLatex(); r.SetNDC(); r.SetTextFont(42); r.SetTextSize(0.04); r.SetTextAlign(31)
    r.DrawLatex(1 - pad.GetRightMargin(), top + 0.012, energy)
    # keep references alive on the pad
    pad._cms_labels = [l, r] + ([l2] if sim else [])


def save(canvas, path_noext, exts=("png", "pdf")):
    for e in exts:
        canvas.SaveAs("%s.%s" % (path_noext, e))
