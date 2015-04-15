from ROOT import *
import os,array,sys
from rebin_tools import *
from slice_tools import *
# from rootpy.interactive import wait

doslices = True
do2d = True
datadir = os.getcwd()
mctypes = ["nl0","nl1","nl2"]
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
  pad.SetTopMargin(0.05)

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
  gev_vars = ["mj", "leadfjm", "leadfjpt", "met", "mt", "isr", "isr1pt", "isr2pt", "isr3pt", "ptt", "leadpttop", "subleadpttop", "leadptglu"]
  if var in gev_vars or var[0:1]=="ht":
    return "GeV"
  else:
    return ""

def set_pallete(samp):
  if (samp=="ttbar"):
    Red = array.array('d',[1., 0.39])
    Green = array.array('d',[1.,0.51])
    Blue = array.array('d',[1., 0.78])
    Length = array.array('d',[0.,1.])
    nb = 50
    TColor.CreateGradientColorTable(2,Length,Red,Green,Blue,nb)
  elif (samp=="T1tttt1500"):
    Red = array.array('d',[1., 0.23])
    Green = array.array('d',[1.,0.72])
    Blue = array.array('d',[1., 0.52])
    Length = array.array('d',[0.,1.])
    nb = 50
    TColor.CreateGradientColorTable(2,Length,Red,Green,Blue,nb)
  elif (samp=="T1tttt1200"):
    Red = array.array('d',[1., 0.78])
    Green = array.array('d',[1.,0.47])
    Blue = array.array('d',[1., 0.63])
    Length = array.array('d',[0.,1.])
    nb = 50
    TColor.CreateGradientColorTable(2,Length,Red,Green,Blue,nb)


gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetPalette(55)
ROOT.gStyle.SetPaintTextFormat(".0f")
gROOT.SetBatch(kTRUE)

# ------------- SELECTIONS -------------
selns = [
  "nl1-ht750-met250-mt150-nj6-nb2",
  # "nl1-ht500-met200-mt150-nj0-nb0",
  # "nl1-ht500-met200-mt150-nj6-nb0",
  # "nl1-ht500-met200-mt0-nj6-nb0",
  # "nlna-ht500-met200-mt0-nj6-nb0",
  # "nlna-ht0-met200-mt0-nj0-nb0",
  # "nlna-ht0-met0to100-mt0-nj0-nb0",
  # "nlna-ht500-met0-mt0-nj0-nb0"
  # "nlna-ht500-met200-mt0-nj6-nb0-mj300to450",
  # "nlna-ht500-met200-mt0-nj6-nb0-mj450to600",
  # "nlna-ht500-met200-mt0-nj6-nb0-mj600toInf"
]

# ------------- SAMPLES -------------
sampldict = {}
sampldict["ttbar"] = "t#bar{t}"
# sampldict["qcd"] = ("QCD",kRed, 3017, "qcd")
# sampldict["wjets"] = ("W+jets", kViolet, 3490, "w")
# sampldict["T1tttt1500"] = "T1tttt NC"
# sampldict["T1tttt1200"] = "T1tttt C"

# colors = [kRed+1, kOrange, kGreen+1, kBlue+1]
colors = [kMagenta+2, kBlue+1, kGreen+1, kOrange-3, kRed+1]

# markers = [20, kBlue+2, kGreen+2, kOrange+7, kBlack]

# --------- VARIABLES -------------
edgelist = {}
edgelist["mj"] = [300., 450., 600.,1500.]
edgelist["fjm1"] = [0., 50.,100.,150.,200.,250.,300.,350.,400.,500.]
edgelist["met"] = [200.,300., 400.,500.,1500.]
edgelist["met_me"] = [200.,300., 400.,500.,1500.]
edgelist["ptt"] = [0.,200., 400.,600.,1000.]
edgelist["avetoppt"] = [0.,200., 400.,600.,1000.]
edgelist["dphi_tt"] = [x/3.*3.14 for x in range(0,3)]
edgelist["dphi_fjm1_fjm2"] = [x/3.*3.14 for x in range(0,3)]
edgelist["ht"] = [500.,750.,1000.,1500., 2000.,4000.]
edgelist["ht_isr_me"] = [0.,100.,300.,500.,1000.,1500.]
edgelist["ht_fsr"] = [1]
edgelist["ht_part"] = [0.,500.,1000.,1500., 2500.]
edgelist["njets"] = [4, 6, 8, 20]
edgelist["npart"] = [0, 4, 6, 8, 20]
edgelist["nisr_me"] = [-0.5,0.5,1.5,2.5,3.5]
edgelist["nfsr"] = [0,1,2,3]
edgelist["leadpttop"] = [0., 100., 200., 300., 400., 600.]

# --------- VARIABLES -------------
varlbl = {}
varlbl["mj"] = "M_{J}"
varlbl["fjm1"] = "m(J_{1})"
varlbl["ptt"] = "p_{T}(t#bar{t})"
varlbl["avetoppt"] = "p^{ave}_{T}(t)"
varlbl["dphi_tt"] = "#Delta#phi(t,#bar{t})"
varlbl["dphi_fjm1_fjm2"] = "#Delta#phi(m_{j1},m_{j2})"
varlbl["ht"] = "H_{T}"
varlbl["ht_isr_me"] = "H_{T} ISR ME"
varlbl["ht_fsr"] = "H_{T} FSR"
varlbl["ht_part"] = "H_{T} part."
varlbl["njets"] = "n_{jets}"
varlbl["npart"] = "n_{partons}"
varlbl["nisr_me"] = "n_{ISR ME}"
varlbl["nfsr"] = "n_{FSR}"
varlbl["leadpttop"] = "p_T^{lead}(t)"
varlbl["met"] = "MET"
varlbl["met_me"] = "MET (ME)"

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
# var_pairs.append(["isr","ptt"])
# var_pairs.append(["nisrjets","mj"])
# var_pairs.append(["dphi_tt","ptt"])
# var_pairs.append(["dphitt","avetoppt"])
# var_pairs.append(["minDphiIsrTop","mj"])
# var_pairs.append(["minDphiIsrTop","leadfjm"])
# var_pairs.append(["dphi_tt","dphi_fjm1_fjm2"])

# var_pairs.append(["npart","njets"])
# var_pairs.append(["nisr_me","mj"])
# var_pairs.append(["fjm1","nisr_me"])
# var_pairs.append(["ht_fsr","mj"])
# var_pairs.append(["nisr_me","ht_isr_me"])

# var_pairs.append(["dphi_tt", "nisr_me"])
# var_pairs.append(["dphi_tt", "ht_isr_me"])
# var_pairs.append(["nfsr","mj"])
# var_pairs.append(["nfsr","dphi_tt"])

# var_pairs.append(["ht_part","ht"])
# var_pairs.append(["ht_isr_me","mj"])
# var_pairs.append(["ht_fsr","mj"])
# var_pairs.append(["ht_part","mj"])

# var_pairs.append(["dphi_tt","ht_isr_me"])

var_pairs.append(["met","met_me"])

# var_pairs.append(["pt_gg","pt_isr"])
# var_pairs.append(["ptt","pt_isr"])

for samp in sampldict.keys():
  set_pallete(samp)
  flist = TFile(os.path.join(datadir,"mj_plots_"+samp+".root"),"READ")
  for seln in selns:
    for mctype in mctypes:

      # can = TCanvas("can"+samp+seln+mctype,"can",1000,1000)
      # for var in varlbl:
      #   pad = make_pad()
      #   pad.Draw()
      #   pad.cd()
      #   pad.SetLogy()

      #   # print "Get histogram: ", '_'.join([seln, var, samp, mctype])
      #   if (mctype!=""): hist = flist.Get('_'.join([seln, var, samp, mctype])).Clone()
      #   else: hist = flist.Get('_'.join([seln, var, samp])).Clone()
      #   hist.SetDirectory(0)
      #   hist.SetFillColor(colors[0])
      #   hist.SetFillStyle(colors[0])
      #   hist.SetLineColor(colors[0])
      #   hist.SetLineWidth(4)
      #   style_axis(hist)
      #   hist.Draw("hist")
        
      #   can.Print(os.path.join(outdir, seln+"_"+var+".pdf"))


      # ---------------------2D distributions -------------------------
      if not do2d: continue
      for pair in var_pairs:
        print samp, seln, mctype, pair
        if (mctype!=""): hist = flist.Get('_'.join([seln, pair[0], pair[1], samp, mctype])).Clone()
        else: hist = flist.Get('_'.join([seln, pair[0], pair[1], samp])).Clone()

        can = TCanvas(hist.GetName(),"can",1000,900)
        pad = make_pad()
        pad.Draw()
        pad.cd()

        slices = {}
        if (doslices):
          slices["X"] = slice2DX(hist, edgelist[pair[0]],pair[0])
          slices["Y"] = slice2DY(hist, edgelist[pair[1]],pair[1])

        style_axis(hist)
        hist.SetLineWidth(3)
        hist.SetLineColor(kRed+1)
        hist.SetMarkerSize(1)
        # hist.RebinX(2)
        if (pair[0]=="mj"): hist.RebinX(2)
        elif (pair[1]=="mj"): hist.RebinY(2)
        if (pair[0]=="ht_isr_me"): hist.RebinX(2)
        elif (pair[1]=="ht_isr_me"): hist.RebinY(2)
        hist.Draw("colz text")
        # hist.Draw("scat")
        # label(0.77,0.85,"#rho={:.2f}".format(hist.GetCorrelationFactor()),1,0.05)
        can.SetLogy()
        can.Print(hist.GetName()+".pdf")

        if (doslices):
          for iaxis,axis in enumerate(["X","Y"]):

            leg = make_legend()
            # leg.SetHeader(seln)

            svar = pair[iaxis]
            ymax = 0
            for i,slice in enumerate(slices[axis]):
              if (ymax < 1.5*slice.GetMaximum()): 
                ymax = 1.5*slice.GetMaximum()
                slices[axis][0].GetYaxis().SetRangeUser(0.,ymax)
              if (i<len(colors)): slice.SetLineColor(colors[i])
              slice.SetLineWidth(3)
              slice.SetLineStyle(len(slices[axis])-i)
              slice.SetYTitle("Events")
              # slice.SetMarkerStyle(20+i)
              # slice.SetMarkerSize(2.5)
              # slice.SetMarkerColor(colors[i])
              # slice.Rebin(2)
              style_axis(slice)
              if (i==0): slice.Draw("hist")
              else: slice.Draw("hist same")
              if (svar == "dphi_tt"):
                leg.AddEntry(slice, '{} {:.1f} - {:.1f} {} (#mu {:.1f}, RMS {:.1f})'.format(varlbl[svar], edgelist[svar][i], edgelist[svar][i+1], units(svar), slice.GetMean(), slice.GetRMS()), "L")
              elif (svar[0]=="n"):
                leg.AddEntry(slice, '{} = {:.0f} {} (#mu {:.1f}, RMS {:.1f})'.format(varlbl[svar], edgelist[svar][i]+0.5, units(svar), slice.GetMean(), slice.GetRMS()), "L")
              else:
                leg.AddEntry(slice, '{} {:.0f} - {:.0f} {} (#mu {:.1f}, RMS {:.1f})'.format(varlbl[svar], edgelist[svar][i], edgelist[svar][i+1], units(svar), slice.GetMean(), slice.GetRMS()), "L")

            leg.Draw()
            
            if (mctype!=""): can.Print(os.path.join(outdir, '_'.join([seln, pair[1-iaxis], "inSlices", svar, samp, mctype+".pdf"])))
            else: can.Print(os.path.join(outdir, '_'.join([seln, pair[1-iaxis], "inSlices", svar, samp+".pdf"])))

            # can.SetLogy(1)
            # if (mctype!=""): can.Print(os.path.join(outdir, '_'.join([seln, pair[1-iaxis], "inSlices", svar, samp, mctype+"_log.pdf"])))
            # else: can.Print(os.path.join(outdir, '_'.join([seln, pair[1-iaxis], "inSlices", svar, samp+"_log.pdf"])))

