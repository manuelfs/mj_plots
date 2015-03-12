from ROOT import TH1F, TH2F, TCanvas
import array
import math

def slice2DX(hist, edge_list, var):
  slices = []
  for i in range(0,len(edge_list)-1):
    low = hist.GetXaxis().FindBin(edge_list[i])
    high = hist.GetXaxis().FindBin(edge_list[i+1]*0.9999) #want the bin that goes up to the edge, not one that starts at edge
    hnm = hist.GetName()
    hnm.replace("_"+var+"_","_")
    slices.append( hist.ProjectionY('_'.join([hnm, var+str(edge_list[i]), str(edge_list[i+1])]), low, high) )
  return slices
  
def slice2DY(hist, edge_list, var):
  slices = []
  for i in range(0,len(edge_list)-1):
    low = hist.GetYaxis().FindBin(edge_list[i])
    high = hist.GetYaxis().FindBin(edge_list[i+1]*0.9999) #want the bin that goes up to the edge, not one that starts at edge
    hnm = hist.GetName()
    hnm.replace("_"+var+"_","_")
    slices.append( hist.ProjectionX('_'.join([hnm, var+str(edge_list[i]), str(edge_list[i+1])]), low, high) )
  return slices

