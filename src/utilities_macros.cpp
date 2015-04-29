//----------------------------------------------------------------------------
// utilities_macros - Various functions used accross the code
//----------------------------------------------------------------------------

#ifndef INT_ROOT
#include "utilities_macros.hpp"
#endif

#include <cmath>
#include <vector>

#include <iostream>
#include <string>
#include <stdexcept>

#include "TString.h"
#include "TH1D.h"
#include "TMath.h"

using namespace std;


hfeats::hfeats(TString ivarname, int inbins, float iminx, float imaxx, vector<int> isamples,
	       TString ititle, TString icuts, float icut, TString itagname):
  title(ititle),
  varname(ivarname),
  cuts(icuts),
  nbins(inbins),
  minx(iminx),
  maxx(imaxx),
  cut(icut),
  samples(isamples),
  tagname(itagname){
  format_tag();
  unit = "";
  if(title.Contains("GeV)")) unit = "GeV";
  if(title.Contains("phi")) unit = "rad";
  }

hfeats::hfeats(TString ivarname, int inbins, float *ibinning, vector<int> isamples,
	       TString ititle, TString icuts, float icut, TString itagname):
  title(ititle),
  varname(ivarname),
  cuts(icuts),
  nbins(inbins),
  binning(ibinning),
  cut(icut),
  samples(isamples),
  tagname(itagname) {
  minx = binning[0]; maxx = binning[nbins];
  format_tag();
  unit = "";
  if(title.Contains("GeV)")) unit = "GeV";
  if(title.Contains("phi")) unit = "rad";
  }

void hfeats::format_tag(){
  tag = varname; 
  if(cuts!="1")   tag+="_"+cuts;
  if(tagname!="") tag+="_"+tagname;

  tag.ReplaceAll(".",""); 
  tag.ReplaceAll("(",""); tag.ReplaceAll("$","");  tag.ReplaceAll(")",""); 
  tag.ReplaceAll("[",""); tag.ReplaceAll("]",""); tag.ReplaceAll("||","_");
  tag.ReplaceAll("/","_"); tag.ReplaceAll("*",""); tag.ReplaceAll("&&","_");
  tag.ReplaceAll(">=","ge"); tag.ReplaceAll("<=","se"); 
  tag.ReplaceAll(">","g"); tag.ReplaceAll("<","s"); tag.ReplaceAll("=","");
  tag.ReplaceAll("+",""); 
}

sfeats::sfeats(vector<TString> ifile, TString ilabel, int icolor, int istyle, TString icut){
  file = ifile; label = ilabel; cut = icut;
  color = icolor; style = istyle;
  isSig = ifile[0].Contains("T1tttt");
  factor = "1";
  tag = label; 
  tag.ReplaceAll("(",""); tag.ReplaceAll(",","_");  tag.ReplaceAll(")",""); 
  tag.ReplaceAll("{",""); tag.ReplaceAll("#,","");  tag.ReplaceAll("}",""); 
}

// Function that calculates the chi2 of a histogram with respect to the flat hypothesis
void calc_chi2(TH1D *histo, float &chi2, int &ndof, float &pvalue, float &average){
  vector<double> vals[2];
  double sumx_sig2(0), sum_sig2(0);
  ndof = 0;
  for(int bin(1); bin<=histo->GetNbinsX(); bin++){
    if(histo->GetBinError(bin) > 0){
      vals[0].push_back(histo->GetBinContent(bin));
      vals[1].push_back(histo->GetBinError(bin));
      
      sumx_sig2 += vals[0][ndof]/pow(vals[1][ndof],2);
      sum_sig2 += 1/pow(vals[1][ndof],2);
      ndof++;
    }
  }
  if(sum_sig2<=0){cout<<"All errors in histo are zero. Exiting."<<endl; return;}
  average = sumx_sig2/sum_sig2;
  chi2 = 0; ndof--;
  for(int ival(0); ival <= ndof; ival++){
    chi2 += pow((vals[0][ival]-average)/vals[1][ival],2);
  }
  pvalue = TMath::Prob(chi2,ndof);
}

// Function that calculates the chi2 from the difference of 2 histos, after normalizing their areas
void calc_chi2_diff(TH1D *histo1, TH1D *histo2, float &chi2, int &ndof, float &pvalue, float *average){
  TH1D *histos[] = {histo1, histo2};
  vector<double> vals[2][2];
  double sumx_sig2[] = {0,0}, sum_sig2[] = {0,0};
  int ndofs[] = {0,0};
  for(int his(0); his<2; his++){
    for(int bin(1); bin<=histos[his]->GetNbinsX(); bin++){
      if(histos[his]->GetBinError(bin) > 0){
	vals[his][0].push_back(histos[his]->GetBinContent(bin));
	vals[his][1].push_back(histos[his]->GetBinError(bin));
	
	sumx_sig2[his] += vals[his][0][ndofs[his]];
	sum_sig2[his] += 1;
	ndofs[his]++;
      }
    }
    if(sum_sig2[his]<=0){cout<<"All errors in histo are zero. Exiting."<<endl; return;}
    average[his] = sumx_sig2[his]/sum_sig2[his];
    ndofs[his]--;
  } // Loop over histos

  // In principle, we should use just bins where both histos have entries. To be modified
  if(ndofs[0] != ndofs[1]) {cout<<"First histo has "<<ndofs[0]<<" ndof and second "<<ndofs[1]<<endl; return;}
  else ndof = ndofs[0];
  chi2 = 0; 
  double Raver = average[0]/average[1];
  for(int ival(0); ival <= ndof; ival++){    
    double error(sqrt(pow(vals[0][1][ival],2)+pow(vals[1][1][ival]*Raver,2)));
    chi2 += pow((vals[0][0][ival]-vals[1][0][ival]*Raver)/error,2);
  }
  pvalue = TMath::Prob(chi2,ndof);
}

