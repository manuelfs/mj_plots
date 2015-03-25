from ROOT import TH1F, TH2F, TCanvas
import array
import math
import sys

def slice2DX(hist, edge_list, var):
  slices = []
  low, high = 0, 0
  for i in range(0,len(edge_list)-1):
    low = hist.GetXaxis().FindBin(edge_list[i]*1.0001)
    print low, hist.GetXaxis().GetBinLowEdge(low+1), edge_list[i], edge_list[i+1]
    if (i>0 and low<=high):
      print "FATAL:: slice_tools: Slices overlap, i.e. edge_list contains values that do not correspond to a bin edge!"
      print low, high
      sys.exit(1)
    high = hist.GetXaxis().FindBin(edge_list[i+1]*0.9999) #want the bin that goes up to the edge, not one that starts at edge
    hnm = hist.GetName()
    hnm.replace("_"+var+"_","_")
    slices.append( hist.ProjectionY('_'.join([hnm, var+str(edge_list[i]), str(edge_list[i+1])]), low, high) )
  return slices
  
def slice2DY(hist, edge_list, var):
  slices = []
  low, high = 0, 0
  for i in range(0,len(edge_list)-1):
    low = hist.GetYaxis().FindBin(edge_list[i])
    if (i>0 and low<=high):
      print "FATAL:: slice_tools: Slices overlap, i.e. edge_list contains values that do not correspond to a bin edge!"
      sys.exit(1)
    high = hist.GetYaxis().FindBin(edge_list[i+1]*0.9999) #want the bin that goes up to the edge, not one that starts at edge
    hnm = hist.GetName()
    hnm.replace("_"+var+"_","_")
    slices.append( hist.ProjectionX('_'.join([hnm, var+str(edge_list[i]), str(edge_list[i+1])]), low, high) )
  return slices

