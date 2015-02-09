from ROOT import TH1F, TH2F, TCanvas
import array
import math

binsLowEdgeListInit = [10., 20.,50.,75., 100., 150., 200., 250.,300., 350.,400.]

#---------------------------------------
#       1D REBIN FUNCITON
#---------------------------------------

def rebin1D(hist, binsLowEdgeList = binsLowEdgeListInit):
  if (len(binsLowEdgeList)==1):
    hnew = hist.Clone(hist.GetName() + "_rebin")
    hnew.Rebin(binsLowEdgeList[0])
  else:
    binsLowEdgeArr = array.array('d',binsLowEdgeList)
    nBinsNew = len(binsLowEdgeList) - 1
    hnew = hist.Rebin(nBinsNew, hist.GetName() + "_rebin", binsLowEdgeArr)

  return hnew

def rebin1D_wOverflow(hist, binsLowEdgeList = binsLowEdgeListInit):
  if (len(binsLowEdgeList)==1):
    hnew = hist.Clone(hist.GetName() + "_rebin")
    hnew.Rebin(binsLowEdgeList[0])
  else:
    binsLowEdgeArr = array.array('d',binsLowEdgeList)
    nBinsNew = len(binsLowEdgeList) - 1
    hnew = hist.Rebin(nBinsNew, hist.GetName() + "_rebin", binsLowEdgeArr)

  lastBin = hnew.GetBinContent(nBinsNew) + hnew.GetBinContent(nBinsNew+1)
  lastBinErr = math.pow(hnew.GetBinError(nBinsNew),2) + math.pow(hnew.GetBinError(nBinsNew+1),2)
  lastBinErr = math.sqrt(lastBinErr)
  hnew.SetBinContent(nBinsNew, lastBin)
  hnew.SetBinError(nBinsNew, lastBinErr)

  return hnew

#---------------------------------------
#       2D REBIN FUNCITON
#---------------------------------------

def rebin2D(hist, binsLowEdgeList = binsLowEdgeListInit):

  # binsLowEdgeList = [300., 350., 400., 450., 500., 550., 650., 750., 1200.]
  binsLowEdgeArr = array.array('d',binsLowEdgeList)
  nBinsNew = len(binsLowEdgeList) - 1

  hnew = TH2F(hist.GetName() + "_rebin", hist.GetTitle() + "_rebin", nBinsNew, binsLowEdgeArr, nBinsNew, binsLowEdgeArr)
  hnew.GetXaxis().SetTitle(hist.GetXaxis().GetTitle())
  hnew.GetYaxis().SetTitle(hist.GetYaxis().GetTitle())
  hnew.Sumw2()

  for i in range(0,nBinsNew):
    binWidthOld = hist.GetBinWidth(1)
    firstBinX = int((binsLowEdgeList[i]-hist.GetXaxis().GetXmin())/binWidthOld) + 1    # 300/50 = 6 +1 . 7
    lastBinX = int((binsLowEdgeList[i+1]-hist.GetXaxis().GetXmin())/binWidthOld)       # 400/50 = 8
    for j in range(0,nBinsNew):
      firstBinY = int((binsLowEdgeList[j]-hist.GetYaxis().GetXmin())/binWidthOld) + 1
      lastBinY = int((binsLowEdgeList[j+1]-hist.GetYaxis().GetXmin())/binWidthOld)
      hnew.SetBinContent(i+1,j+1, hist.Integral(firstBinX, lastBinX, firstBinY, lastBinY))

      err2 = 0.
      for iBinX in range(firstBinX, lastBinX+1):
        for iBinY in range(firstBinY, lastBinY+1):
          err2 += pow(hist.GetBinError(iBinX,iBinY),2)
        
      hnew.SetBinError(i+1,j+1, math.sqrt(err2))

  return hnew

#---------------------------------------
#       2D REBIN FUNCITON
#---------------------------------------

def rebin2D_Xonly(hist, binsLowEdgeList = binsLowEdgeListInit):

  binsLowEdgeArr = array.array('d',binsLowEdgeList)
  nBinsNew = len(binsLowEdgeList) - 1

  highedge = hist.GetYaxis().GetBinLowEdge(hist.GetNbinsY()) + hist.GetYaxis().GetBinWidth(hist.GetNbinsY())  
  hnew = TH2F(hist.GetName() + "_rebin", hist.GetTitle() + "_rebin", nBinsNew, binsLowEdgeArr, hist.GetNbinsY(), hist.GetYaxis().GetBinLowEdge(1),highedge)
  hnew.GetXaxis().SetTitle(hist.GetXaxis().GetTitle())
  hnew.GetYaxis().SetTitle(hist.GetYaxis().GetTitle())
  hnew.Sumw2()

  for i in range(0,nBinsNew):
    binWidthOld = hist.GetBinWidth(1)
    firstBinX = int(binsLowEdgeList[i]/binWidthOld) + 1    # 300/50 = 6 +1 . 7
    lastBinX = int(binsLowEdgeList[i+1]/binWidthOld)       # 400/50 = 8
    for j in range(0,hist.GetNbinsY()):
      hnew.SetBinContent(i+1,j+1, hist.Integral(firstBinX, lastBinX, j+1, j+1))

      err2 = 0.
      for iBinX in range(firstBinX, lastBinX+1):
        err2 += pow(hist.GetBinError(iBinX,j+1),2)
        
      hnew.SetBinError(i+1,j+1, math.sqrt(err2))

  return hnew

# Function that makes bin heights independent of bin widths
def smoothBinning(hist):
  for i in range(0,hist.GetNbinsX()):
    denom = float(hist.GetBinWidth(i+1))
    newBin = float(hist.GetBinContent(i+1))/denom
    newErr = float(hist.GetBinError(i+1))/denom
    hist.SetBinContent(i+1,newBin)
    hist.SetBinError(i+1,newErr)
  
  return
