from ROOT import *
import os,array,sys
from rebin_tools import *
from slice_tools import *
# from rootpy.interactive import wait

datadir = os.getcwd()
print "Looking for input in", datadir

outdir = os.getcwd()
if (not os.path.exists(outdir)):
  os.mkdir(outdir)
print "Output directory set to",outdir

def make_pad():
  pad = TPad("pad"," ",0.,0.,1.,1.)
  pad.SetFillColor(0)
  pad.SetLeftMargin(0.17)
  pad.SetRightMargin(0.05)
  pad.SetBottomMargin(0.15)
  pad.SetTopMargin(0.02)

  if (type(hist)==type(TH2D())): pad.SetRightMargin(0.1)

  return pad

def make_legend():
  legX, legY = 0.2, 0.95
  legW, legH = 0.3, 0.07*3.
  leg = TLegend(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(0.03); 
  leg.SetTextSize(0.04); 
  leg.SetFillColor(0); 
  leg.SetFillStyle(0); 
  leg.SetBorderSize(0);
  leg.SetTextFont(42);
  return leg

def style_axis(hist):
  ttl_size = 0.06
  ttl_offset = 1.3

  if (type(hist)!=type(TH3D())): ttl_offset = ttl_offset - 0.2

  hist.GetXaxis().SetLabelSize(ttl_size - 0.01)
  hist.GetXaxis().SetTitleSize(ttl_size)
  hist.GetXaxis().SetTitleOffset(ttl_offset-0.2)
  hist.GetYaxis().SetLabelSize(ttl_size - 0.01)
  hist.GetYaxis().SetTitleSize(ttl_size)
  hist.GetYaxis().SetTitleOffset(ttl_offset+0.2)
  hist.GetZaxis().SetLabelSize(ttl_size - 0.01)
  hist.GetZaxis().SetTitleSize(ttl_size)
  hist.GetZaxis().SetTitleOffset(ttl_offset)


def label(x, y, text, color=1, size=0.08): 
  l=TLatex()
  l.SetTextSize(size); 
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);

def units(var): 
  gev_vars = ["mj", "ht", "leadfjm", "leadfjpt", "met", "mt", "isr", "isr1pt", "isr2pt", "isr3pt", "ptt", "leadpttop", "subleadpttop", "leadptglu"]
  if var in gev_vars:
    return "GeV"
  else:
    return ""


gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetPalette(55)
ROOT.gStyle.SetPaintTextFormat(".1f");

Red = array.array('d',[1., 0.39])
Green = array.array('d',[1.,0.51])
Blue = array.array('d',[1., 0.78])
Length = array.array('d',[0.,1.])
nb = 50
TColor.CreateGradientColorTable(2,Length,Red,Green,Blue,nb)

gROOT.SetBatch(kTRUE)
# ROOT.gStyle.SetPaintTextFormat(".0f");

# ------------- SELECTIONS -------------
selns = [
  "nl1-ht750-met250-mt150-nj6-nb2",
  "nl1-ht500-met200-mt150-nj0-nb0",
  "nl1-ht500-met200-mt150-nj6-nb0",
  "nl1-ht500-met200-mt0-nj6-nb0",
  "nl1-ht500-met200-mt0-nj6-nb0-mj0to200",
  "nl1-ht500-met200-mt0-nj6-nb0-mj200to400",
  "nl1-ht500-met200-mt0-nj6-nb0-mj400to600",
  "nl1-ht500-met200-mt0-nj6-nb0-mj600toInf"
]

# ------------- SAMPLES -------------
sampldict = {}
sampldict["ttbar"] = "t#bar{t}"
# sampldict["qcd"] = ("QCD",kRed, 3017, "qcd")
# sampldict["wjets"] = ("W+jets", kViolet, 3490, "w")
# sampldict["T1tttt1500"] = ("T1tttt NC", kGreen+2, 3013,"t1tn")
# sampldict["T1tttt1200"] = ("T1tttt C", kBlack, 3013,"t1tc")

colors = [kRed+1, kBlue+2, kGreen+2, kOrange+7, kBlack]

# ------ Samples used in plots depend on what is defined in sampldict (not samp_order)
samp = 'ttbar'

# --------- VARIABLES -------------
edgelist = {}
edgelist["mj"] = [0.,200.,300., 400.,600.,1500.]
edgelist["ptt"] = [0.,200., 400.,600.,1000.]
edgelist["avetoppt"] = [0.,200., 400.,600.,1000.]
edgelist["dphi_tt"] = [x/5.*3.14 for x in range(0,5)]
edgelist["dphi_fjm1_fjm2"] = [x/5.*3.14 for x in range(0,5)]
edgelist["ht"] = [500.,750.,1000.,1500., 2000.,4000.]
edgelist["njets"] = [4, 6, 8, 20]
edgelist["nisrjets"] = [0,1,2,3,4]
edgelist["leadpttop"] = [0., 100., 200., 300., 400., 600.]

# --------- VARIABLES -------------
varlbl = {}
varlbl["mj"] = "M_{J}"
varlbl["ptt"] = "p_{T}(t#bar{t})"
varlbl["avetoppt"] = "p^{ave}_{T}(t)"
varlbl["dphi_tt"] = "#Delta#phi(t,#bar{t})"
varlbl["dphi_fjm1_fjm2"] = "#Delta#phi(m_{j1},m_{j2})"
varlbl["ht"] = "H_{T}"
varlbl["njets"] = "n_{jets}"
varlbl["nisrjets"] = "n_{ISR}"
varlbl["leadpttop"] = "p_T^{lead}(t)"

# --------- PAIRS OF VARIABLES TO PLOT -------------
var_pairs = []
# var_pairs.append(["ht","mj"])
# var_pairs.append(["njets","mj"])
# var_pairs.append(["njets","ht"])
# var_pairs.append(["nbl","mj"])
# var_pairs.append(["met","mj"])
# var_pairs.append(["met","ht"])
# var_pairs.append(["met","njets"])
# var_pairs.append(["ptt","ht"])
# var_pairs.append(["leadpttop","ht"])

# var_pairs.append(["leadpttop","ptt"])
# var_pairs.append(["ptt","mj"])
# var_pairs.append(["leadpttop","mj"])
# var_pairs.append(["dphitt","mj"])
# var_pairs.append(["leadpttop","dphitt"])
# var_pairs.append(["leadfjm","dphitt"])
# var_pairs.append(["dphitt","nisrjets"])
# var_pairs.append(["isr","ptt"])
# var_pairs.append(["nisrjets","mj"])
# var_pairs.append(["dphitt","ptt"])
# var_pairs.append(["dphitt","avetoppt"])
# var_pairs.append(["isr1pt","avetoppt"])
# var_pairs.append(["isr2pt","avetoppt"])
# var_pairs.append(["isr3pt","avetoppt"])
# var_pairs.append(["minDphiIsrTop","mj"])
# var_pairs.append(["minDphiIsrTop","leadfjm"])
var_pairs.append(["dphi_tt","dphi_fjm1_fjm2"])

flist = TFile(os.path.join(datadir,"mj_plots_"+samp+".root"),"READ")

for seln in selns:
  for pair in var_pairs:
    hist = flist.Get('_'.join([seln, pair[0], pair[1], samp])).Clone()
    hist.RebinY(2)
    slices = {}
    slices["X"] = slice2DX(hist, edgelist[pair[0]],pair[0])
    slices["Y"] = slice2DY(hist, edgelist[pair[1]],pair[1])

    can = TCanvas(hist.GetName(),"can",1000,900)
    pad = make_pad()
    pad.Draw()
    pad.cd()

    hist.SetLineWidth(3)
    hist.SetLineColor(kRed+1)
    hist.SetMarkerSize(1.8)
    style_axis(hist)
    hist.Draw("colz text")
    # label(0.77,0.85,"#rho={:.2f}".format(hist.GetCorrelationFactor()),1,0.05)
    can.Print(hist.GetName()+".pdf")

    for iaxis,axis in enumerate(["X","Y"]):

      leg = make_legend()
      # leg.SetHeader(seln)

      ymax = 0
      for i,slice in enumerate(slices[axis]):
        if (ymax < slice.GetMaximum()): 
          ymax = 1.5*slice.GetMaximum()
          slices[axis][0].GetYaxis().SetRangeUser(0.,ymax)
        slice.SetLineColor(colors[i])
        slice.SetLineWidth(4)
        slice.SetLineStyle(i+2)
        # slice.Rebin(2)
        style_axis(slice)
        if (i==0): slice.Draw("hist")
        else: slice.Draw("hist same")
        svar = pair[iaxis]
        if (svar == "dphitt"):
          leg.AddEntry(slice, '{} {:.1f} - {:.1f} {} (#mu {:.1f}, RMS {:.1f})'.format(varlbl[svar], edgelist[svar][i], edgelist[svar][i+1], units(svar), slice.GetMean(), slice.GetRMS()), "F")
        else:
          leg.AddEntry(slice, '{} {:.0f} - {:.0f} {} (#mu {:.1f}, RMS {:.1f})'.format(varlbl[svar], edgelist[svar][i], edgelist[svar][i+1], units(svar), slice.GetMean(), slice.GetRMS()), "F")

      leg.Draw()
      
      can.Print(os.path.join(outdir, seln + "_"+pair[1-iaxis]+"_inSlices_"+svar+"_"+samp+".pdf"))


