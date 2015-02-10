// skim_ntuples.C: Skims reduced trees

#define INT_ROOT
#include "inc/utilities.hpp"
#include "src/utilities.cpp"

#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "TSystem.h"
#include "TDirectory.h"

using namespace std;
using std::cout;
using std::endl;

void skim_ntuples(TString infiles, TString cuts="ht>=500&&met>=200"){
  TString folder(infiles), outfile(infiles);
  folder.Remove(folder.Last('/')+1, folder.Length());
  TString outfolder=folder+"skim/";

  // Finding outfile name
  outfile.Remove(0, outfile.Last('/')); outfile.ReplaceAll("*","");
  if(outfile.Contains(".root")) outfile.ReplaceAll(".root","_"+cuts+".root");
  else outfile += ("_"+cuts+".root");
  outfile.ReplaceAll(">=","ge"); outfile.ReplaceAll("<=","se"); outfile.ReplaceAll("&&","_");
  outfile.ReplaceAll(">","g"); outfile.ReplaceAll("<","s"); outfile.ReplaceAll("=","");
  outfile = outfolder+outfile;

  gSystem->mkdir(outfolder, kTRUE);
  TFile out_rootfile(outfile, "RECREATE");
  out_rootfile.cd();
  TChain tree("tree"); 
  int nfiles = tree.Add(infiles); 

  cout<<"Skimming the "<<nfiles<<" files in "<<infiles<<endl;
  TTree *ctree = tree.CopyTree(cuts);
  ctree->Write();
  out_rootfile.Close();
  cout<<"Written "<<outfile<<endl;

}

void all_skim(TString folder){
  vector<TString> files = dirlist(folder, ".root");

  for(unsigned file(0); file < files.size(); file++){
    skim_ntuples(folder+files[file], "ht30>=500&&met>=200");
  }
}
