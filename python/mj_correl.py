from ROOT import *
import os,array,sys
from rebin_tools import *
# from rootpy.interactive import wait

mode30 = False
# mode30 = True
datadir = os.path.join(os.getcwd(),"out/jet40/")
if (mode30): datadir = "./out/jet30/"

outdir = os.path.join(os.getcwd(),"out/jet40/nfjets/")
if (not os.path.exists(outdir)):
  os.mkdir(outdir)

def make_pad():
  pad = TPad("pad"," ",0.,0.,1.,1.)
  pad.SetFillColor(0)
  pad.SetLeftMargin(0.15)
  pad.SetRightMargin(0.05)
  pad.SetBottomMargin(0.15)
  return pad

def make_legend():
  legX, legY = 0.5, 0.88
  legW, legH = 0.3, 0.07*3.
  leg = TLegend(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(0.048); 
  leg.SetFillColor(0); 
  leg.SetFillStyle(0); 
  leg.SetBorderSize(0);
  leg.SetTextFont(42);
  return leg

def style_axis(hist):
  ttl_size = 0.05
  ttl_offset = 1.3
  if (mode30): 
    ttl_size = 0.04
    ttl_offset = 1.05

  if (type(hist)!=type(TH3D())): ttl_offset = ttl_offset - 0.2

  hist.GetXaxis().SetLabelSize(0.04)
  hist.GetXaxis().SetTitleSize(ttl_size)
  hist.GetXaxis().SetTitleOffset(ttl_offset)
  hist.GetYaxis().SetLabelSize(0.04)
  hist.GetYaxis().SetTitleSize(ttl_size)
  hist.GetYaxis().SetTitleOffset(ttl_offset+0.2)
  hist.GetZaxis().SetLabelSize(0.04)
  hist.GetZaxis().SetTitleSize(ttl_size)
  hist.GetZaxis().SetTitleOffset(ttl_offset)


gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gROOT.SetBatch(kTRUE)
# ROOT.gStyle.SetPaintTextFormat(".0f");

# ------------- SELECTIONS -------------
selndict = {"2L-base": "2-lepton", "1L-base": "1-lepton", "0L-base": "0-lepton"}

# ------------- SAMPLES -------------
sampldict = {}
sampldict["ttbar"] = ("t#bar{t}", kRed+1,3008)
# sampldict["qcd"] = ("QCD",kOrange, 3017)
sampldict["wjets"] = ("W+jets", kBlue+1, 3490)
sampldict["T1tttt1500"] = ("T1t.NC", kGreen+1, 3013)
# sampldict["T1tttt1200"] = ("T1t.C", kBlack, 3013)
# sampldict["T1qqqq1400"] = ("T1q.NC", kBlue, 3013)
# sampldict["T1qqqq1000"] = ("T1q.C", kOrange, 3013)

# ------ Samples used in plots depend on what is defined in sampldict (not samp_order)
samp_order = []
samp_order.append("wjets")
samp_order.append("ttbar")
samp_order.append("qcd")
samp_order.append("T1tttt1500")
samp_order.append("T1tttt1200")
samp_order.append("T1qqqq1400")
samp_order.append("T1qqqq1000")

# ------ VARIABLES TO USE IN THE PLOTS -------------
varlist = []
varlist.append("met")
varlist.append("mt")
varlist.append("mindphin_metjet")
if (mode30):
  varlist.append("mj_30")
  varlist.append("ht30")
  varlist.append("njets30")
  varlist.append("nbl30")
else:
  varlist.append("mj_40")
  varlist.append("mj_40_fjm70")
  varlist.append("ht")
  varlist.append("njets")
  varlist.append("nfjets_40")
  varlist.append("fjets_40_m")
  varlist.append("nbl")

# --------- PAIRS OF VARIABLES TO PLOT -------------
var_pairs = []
if (mode30):
  var_pairs.append(["ht30","mj_30"])
  var_pairs.append(["njets30","mj_30"])
  var_pairs.append(["nbl30","mj_30"])
  var_pairs.append(["met","mj_30"])
  var_pairs.append(["met","njets30"])
  var_pairs.append(["mt","mj_30"])
  var_pairs.append(["mindphin_metjet","mj_30"])
else:
  var_pairs.append(["ht","mj_40"])
  var_pairs.append(["njets","mj_40"])
  var_pairs.append(["njets","ht"])
  var_pairs.append(["nbl","mj_40"])
  var_pairs.append(["met","mj_40"])
  var_pairs.append(["met","ht"])
  var_pairs.append(["met","njets"])
  var_pairs.append(["mt","mj_40"])
  var_pairs.append(["mindphin_metjet","mj_40"])

# --------- SETS OF VARIABLES TO PLOT -------------
var_sets = []
if (mode30):
  var_sets.append(["mj_30","njets30","met"])
  var_sets.append(["ht30","njets30","met"])
else:
  var_sets.append(["mj_40","njets","met"])
  var_sets.append(["ht","njets","met"])

flist = {}
for samp in sampldict.keys():
  flist[samp] = TFile(datadir+"mj_plots_"+samp+".root","READ")

for seln in selndict.keys():
  #------- 1D - DISTRIBUTIONS ---------
  for var in varlist:
    leg = make_legend()
    leg.SetHeader(selndict[seln]+" *SHAPE*")

    can = TCanvas("can"+seln+var,"can",1000,1000)
    pad = make_pad()
    pad.Draw()
    pad.cd()

    hist = {}
    first_draw = True
    for i,samp in enumerate(samp_order):
      if samp not in sampldict.keys(): continue
      hist[samp] = flist[samp].Get(seln+"_"+var+"_"+samp).Clone()
      hist[samp].SetDirectory(0)
      hist[samp].Scale(1./hist[samp].Integral())

      hist[samp].SetFillColor(sampldict[samp][1])
      hist[samp].SetFillStyle(sampldict[samp][2])
      hist[samp].SetLineColor(sampldict[samp][1])
      hist[samp].SetLineWidth(2)

      if (first_draw):
        style_axis(hist[samp])
        hist[samp].GetYaxis().SetRangeUser(0.,0.5)
        hist[samp].Draw("hist")
        first_draw = False
      else: 
        hist[samp].Draw("hist same")

      
    for samp in samp_order:
      if samp not in sampldict.keys(): continue
      leg.AddEntry(hist[samp], sampldict[samp][0], "LP")

    leg.Draw()
    can.Print(outdir+seln+"_"+var+".pdf")

  #------- CORRELATIONS ---------
  for pair in var_pairs:
    if pair[0] not in varlist: continue     
    if pair[1] not in varlist: continue     

    leg = make_legend()
    leg.SetHeader(selndict[seln]+" *SHAPE*")

    can = TCanvas("can"+seln+pair[0]+pair[1],"can",1000,1000)
    pad = make_pad()
    pad.Draw()
    pad.cd()

    hist = {}
    first_draw = True
    for i,samp in enumerate(samp_order):
      if samp not in sampldict.keys(): continue
      if ("wjets" in samp): continue
      hist[samp] = flist[samp].Get(seln+"_"+pair[0]+"_"+pair[1]+"_"+samp).Clone()
      hist[samp].SetDirectory(0)
      hist[samp].Scale(1./hist[samp].Integral())

      hist[samp].SetLineColor(sampldict[samp][1])
      hist[samp].SetLineWidth(2)
      if (first_draw):
        style_axis(hist[samp])
        hist[samp].Draw("box")
        first_draw = False
      else: 
        hist[samp].Draw("box same")

      
    for samp in samp_order:
      if samp not in sampldict.keys(): continue
      if ("wjets" in samp): continue
      corr = hist[samp].GetCorrelationFactor()
      leg.AddEntry(hist[samp], sampldict[samp][0]+" #rho="+"{0:.2f}".format(corr), "F")

    leg.Draw()
    can.Print(outdir+seln+"_"+pair[0]+"_"+pair[1]+".pdf")

  #------- 3D - DISTRIBUTIONS ---------
  for set in var_sets:

    can = TCanvas("can"+seln+set[0]+set[1]+set[2],"can",1000,1000)
    pad = make_pad()
    pad.Draw()
    pad.cd()

    hist = {}
    first_draw = True
    for i,samp in enumerate(samp_order):
      if samp not in sampldict.keys(): continue
      if ("wjets" in samp): continue
      hist[samp] = flist[samp].Get(seln+"_"+set[0]+"_"+set[1]+"_"+set[2]+"_"+samp).Clone()
      hist[samp].RebinX(2)
      hist[samp].RebinY(2)
      hist[samp].RebinZ(2)
      hist[samp].SetDirectory(0)
      hist[samp].Scale(1./hist[samp].Integral())

      hist[samp].SetMarkerColor(sampldict[samp][1])
      if (first_draw):
        style_axis(hist[samp])
        hist[samp].Draw("box")
        first_draw = False
      else: 
        hist[samp].Draw("box same")

    can.Print(outdir+seln+"_"+set[0]+"_"+set[1]+"_"+set[2]+".pdf")



