from ROOT import *
import os,array,sys
from rebin_tools import *
# from rootpy.interactive import wait

mode30 = False
make1d = True
make2d = False
make3d = False

shape = True
# opt = "sign_shape"
# opt = "sign"
# opt = "base_shape"
opt = "mt_cut"
# opt = "mj_ge300"
datadir = os.path.join(os.getcwd(),"out",opt)
print "Looking for input in", datadir

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
selndict = {
  "2L-base": "2l, Preseln", 
  "1L-base": "1l, Preseln", 
  "0L-base": "0l, Preseln",
  "2L-sign": "2l, Preseln,6j,2b", 
  "1L-sign": "1l, Preseln,6j,2b", 
  "0L-sign": "0l, Preseln,6j,2b"
}

# ------------- SAMPLES -------------
sampldict = {}
sampldict["ttbar"] = ("t#bar{t}", kRed+1,3008, "tt")
# sampldict["qcd"] = ("QCD",kRed, 3017, "qcd")
# sampldict["wjets"] = ("W+jets", kViolet, 3490, "w")
sampldict["T1tttt1500"] = ("T1tttt NC", kGreen+2, 3013,"t1tn")
# sampldict["T1tttt1200"] = ("T1tttt C", kBlack, 3013,"t1tc")
# sampldict["T1qqqq1400"] = ("T1qqqq NC", kBlue+1, 3013,"t1qn")
# sampldict["T1qqqq1000"] = ("T1q.C", kOrange, 3013,"t1qc")
# sampldict["T2tt650"] = ("T2t.C", kBlue, 3013,"t2tc") # no skim available
# sampldict["T2tt850"] = ("T2t.NC", kMagenta+3, 3013,")

samplstr = ""
for i,samp in enumerate(sampldict.keys()):
  samplstr = samplstr + sampldict[samp][3]
  if (i<len(sampldict)-1): samplstr = samplstr + "_"

outdir = os.path.join(os.getcwd(),"out",opt,samplstr)
if (not os.path.exists(outdir)):
  os.mkdir(outdir)

# ------ Samples used in plots depend on what is defined in sampldict (not samp_order)
samp_order = []
samp_order.append("wjets")
samp_order.append("ttbar")
samp_order.append("qcd")
samp_order.append("T1tttt1200")
samp_order.append("T1tttt1500")
samp_order.append("T2tt850")
# samp_order.append("T2tt650")
samp_order.append("T1qqqq1400")
samp_order.append("T1qqqq1000")

# ------ VARIABLES TO USE IN THE PLOTS -------------
varlist = []
varlist.append("met")
varlist.append("mt")
varlist.append("mindphin_metjet")
varlist.append("pt_tt")
varlist.append("lead_pt_top")
varlist.append("sublead_pt_top")
varlist.append("lead_pt_gluon")
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
  varlist.append("lead_fjets_40_m")
  varlist.append("lead_fjets_40_pt")
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
else:
  var_pairs.append(["ht","mj_40"])
  var_pairs.append(["njets","mj_40"])
  var_pairs.append(["njets","ht"])
  var_pairs.append(["nbl","mj_40"])
  var_pairs.append(["met","mj_40"])
  var_pairs.append(["met","ht"])
  var_pairs.append(["met","njets"])
  var_pairs.append(["pt_tt","ht"])
  var_pairs.append(["lead_pt_top","ht"])
  var_pairs.append(["pt_tt","mj_40"])
  var_pairs.append(["lead_pt_top","mj_40"])
  var_pairs.append(["lead_pt_top", "pt_tt"])
  var_pairs.append(["isr", "ptt"])

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
  flist[samp] = TFile(os.path.join(datadir,"mj_plots_"+samp+".root"),"READ")

for seln in selndict.keys():
  #------- 1D - DISTRIBUTIONS ---------
  if (make1d):
    for var in varlist:
      leg = make_legend()
      if (shape): leg.SetHeader(selndict[seln]+" *SHAPE*")
      else: leg.SetHeader(selndict[seln])

      can = TCanvas("can"+seln+var,"can",1000,1000)
      pad = make_pad()
      pad.Draw()
      pad.cd()
      if (not shape): pad.SetLogy()

      hist = {}
      first_draw = True
      for i,samp in enumerate(samp_order):
        if samp not in sampldict.keys(): continue
        try:
          hist[samp] = flist[samp].Get(seln+"_"+var+"_"+samp).Clone()
        except:
          print "\n------ Skipping missing:", " seln=", seln, " var=", var, " samp=", samp
          continue
        hist[samp].SetDirectory(0)
        if (shape):
          if (hist[samp].Integral()>0.):
            hist[samp].Scale(1./hist[samp].Integral())

        hist[samp].SetFillColor(sampldict[samp][1])
        hist[samp].SetFillStyle(sampldict[samp][2])
        hist[samp].SetLineColor(sampldict[samp][1])
        hist[samp].SetLineWidth(2)

        if (first_draw):
          style_axis(hist[samp])
          if (shape): 
            hist[samp].GetYaxis().SetRangeUser(0.,0.5)
            hist[samp].GetYaxis().SetTitle("Fraction")
          hist[samp].Draw("hist")
          first_draw = False
        else: 
          hist[samp].Draw("hist same")

      
      if (len(hist)):  
        for samp in samp_order:
          if samp not in hist.keys(): continue
          leg.AddEntry(hist[samp], sampldict[samp][0], "F")

        leg.Draw()
        can.Print(os.path.join(outdir, seln+"_"+var+".pdf"))

  #------- CORRELATIONS ---------
  if (make2d):
    for pair in var_pairs:

      can = TCanvas("can"+seln+pair[0]+pair[1],"can",1000,1000)
      pad = make_pad()
      pad.Draw()
      pad.cd()
      if (not shape): pad.SetLogz()

      hist = {}
      first_draw = True
      for i,samp in enumerate(samp_order):
        if samp not in sampldict.keys(): continue
        try:
          hist[samp] = flist[samp].Get(seln+"_"+pair[0]+"_"+pair[1]+"_"+samp).Clone()
        except:
          print "\n------ Skipping missing:", " seln=", seln, " vars=", pair[0],"-",pair[1], " samp=", samp
          continue
        hist[samp].SetDirectory(0)
        if (hist[samp].Integral()>0.):
          hist[samp].Scale(1./hist[samp].Integral())

        hist[samp].SetLineColor(sampldict[samp][1])
        hist[samp].SetLineWidth(2)
        if (first_draw):
          style_axis(hist[samp])
          # hist[samp].GetZaxis().SetRangeUser(0.,0.)
          hist[samp].Draw("box")
          first_draw = False
        else: 
          hist[samp].Draw("box same")

        
      if (len(hist)):
        for samp in samp_order:
          if samp not in hist.keys(): continue
          corr = hist[samp].GetCorrelationFactor()
          leg.AddEntry(hist[samp], sampldict[samp][0]+" #rho="+"{0:.2f}".format(corr), "F")

        leg.Draw()
        can.Print(os.path.join(outdir, seln+"_"+pair[0]+"_"+pair[1]+".pdf"))

  #------- 3D - DISTRIBUTIONS ---------
  if (make3d):
    for set in var_sets:

      can = TCanvas("can"+seln+set[0]+set[1]+set[2],"can",1000,1000)
      pad = make_pad()
      pad.Draw()
      pad.cd()

      hist = {}
      first_draw = True
      for i,samp in enumerate(samp_order):
        if samp not in sampldict.keys(): continue
        try:
          hist[samp] = flist[samp].Get(seln+"_"+set[0]+"_"+set[1]+"_"+set[2]+"_"+samp).Clone()
        except:
          print "\n------ Skipping missing:", " seln=", seln, " vars=", set[0],"-",set[1],"-",set[2], " samp=", samp
          continue
        hist[samp].RebinX(2)
        hist[samp].RebinY(2)
        hist[samp].RebinZ(2)
        hist[samp].SetDirectory(0)
        if (shape):
          if (hist[samp].Integral()>0.):
            hist[samp].Scale(1./hist[samp].Integral())

        hist[samp].SetMarkerColor(sampldict[samp][1])
        if (first_draw):
          style_axis(hist[samp])
          hist[samp].Draw("box")
          first_draw = False
        else: 
          hist[samp].Draw("box same")

      if (len(hist)):
        can.Print(os.path.join(outdir, seln+"_"+set[0]+"_"+set[1]+"_"+set[2]+".pdf"))



