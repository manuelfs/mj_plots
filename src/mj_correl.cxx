#include "mj_correl.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <set>
#include <map>
#include <vector>

#include "TFile.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLorentzVector.h"

#include "RooStats/NumberCountingUtils.h"

using namespace std;

namespace{
  const double dblmax = numeric_limits<double>::max();
  const int intmax = numeric_limits<int>::max();
}

bool msgdbg = true;
bool mode30 = false;

int main(){

  TString basedir = "/afs/cern.ch/user/m/manuelf/work/ucsb/15-01-30_skim/";
  if (mode30) basedir = "/afs/cern.ch/user/m/manuelf/work/ucsb/15-01-30_skim/ht30/";

  vector<sample> samples;
  samples.push_back(sample("*_TTJets*", "ttbar", true));
  samples.push_back(sample("*-T1tttt_2J_mGl-1500_mLSP-100_*PU20*", "T1tttt1500", true));
  samples.push_back(sample("*-T1tttt_2J_mGl-1200_mLSP-800_*PU20*", "T1tttt1200", true));
  samples.push_back(sample("*-T1qqqq_2J_mGl-1000_mLSP-800_*", "T1qqqq1000", false));
  samples.push_back(sample("*-T1qqqq_2J_mGl-1400_mLSP-100_*", "T1qqqq1400", false));
  samples.push_back(sample("*-T2tt_2J_mStop-650_mLSP-325_*", "T2tt650", true));
  samples.push_back(sample("*-T2tt_2J_mStop-850_mLSP-100_*", "T2tt850", true));
  samples.push_back(sample("*_WJetsToLNu*", "wjets", false));
  samples.push_back(sample("*_QCD_*", "qcd", false));

  //---------------- DEFINE EVENT SELECTIONS ----------------
  // double els_iso_max = 0.1;
  // double mus_iso_max = 0.4;
  
  vector<seln> selns;
  selns.push_back(seln("0L-base",
                /*nleps_*/ 0, /*nisotks_max_*/ intmax,
                /*sig_lep_pt_min*/ 15.,/*veto_lep_pt_min*/ 15.,
                /*ht_min_*/ 500., /*ht_max_*/ dblmax,
                /*mj_min_*/ 0., /*mj_max_*/ dblmax,
                /*met_min_*/ 200., /*met_max_*/ dblmax,
                /*mt_min_*/ 0, /*mt_max_*/ dblmax,
                /*njets_min_*/ 3, /*njets_max_*/ intmax,
                /*nbl_min_*/ 0, /*nbl_max_*/ intmax));
  selns.push_back(seln("1L-base",
                  /*nleps_*/ 1, /*nisotks_max_*/ intmax,
                  /*sig_lep_pt_min*/ 15.,/*veto_lep_pt_min*/ 15.,
                  /*ht_min_*/ 500., /*ht_max_*/ dblmax,
                  /*mj_min_*/ 0., /*mj_max_*/ dblmax,
                  /*met_min_*/ 200., /*met_max_*/ dblmax,
                  /*mt_min_*/ 0, /*mt_max_*/ dblmax,
                  /*njets_min_*/ 3, /*njets_max_*/ intmax,
                  /*nbl_min_*/ 0, /*nbl_max_*/ intmax));

  selns.push_back(seln("2L-base",
                /*nleps_*/ 2, /*nisotks_max_*/ intmax,
                /*sig_lep_pt_min*/ 15.,/*veto_lep_pt_min*/ 15.,
                /*ht_min_*/ 500., /*ht_max_*/ dblmax,
                /*mj_min_*/ 0., /*mj_max_*/ dblmax,
                /*met_min_*/ 200., /*met_max_*/ dblmax,
                /*mt_min_*/ 0, /*mt_max_*/ dblmax,
                /*njets_min_*/ 3, /*njets_max_*/ intmax,
                /*nbl_min_*/ 0, /*nbl_max_*/ intmax));

  //---------------- DEFINE VARIABLES FOR ALL PLOTS ----------------
  map<TString, var> vars_map;
  if (mode30){
    vars_map["mj_30"] = var("mj_30", "Sum(m^{30}_{J}) [GeV]", 30, 0., 1500.);
    vars_map["ht30"] = var("ht30", "H^{30}_{T} [GeV]", 35, 500., 4000.);
    vars_map["njets30"] = var("njets30", "jet^{30} multiplicity", 20, 0., 20.);
    vars_map["nbl30"] = var("nbl30", "b-jet^{30} multiplicity", 6, 0., 6.);
  } else {
    vars_map["mj_40"] = var("mj_40", "Sum(m_{J}) [GeV]", 30, 0., 1500.);
    vars_map["mj_40_fjm70"] = var("mj_40_fjm70", "Sum(m_{J}) (w/ m_{J}>70) [GeV]", 30, 0., 1500.);
    vars_map["ht"] = var("ht", "H_{T} [GeV]", 35, 500., 4000.);
    vars_map["njets"] = var("njets", "jet multiplicity", 20, 0., 20.);
    vars_map["nbl"] = var("nbl", "b-jet multiplicity", 6, 0., 6.);
    vars_map["fjets_40_m"] = var("fjets_40_m", "large-R jet mass [GeV]", 40, 0., 400.);
    vars_map["nfjets_40"] = var("nfjets_40", "large-R jet multiplicity", 20, 0., 20.);
  }
  vars_map["met"] = var("met", "MET [GeV]", 30, 0., 1500.);
  vars_map["mt"] = var("mt", "m_{T} [GeV]", 20, 0., 1000.);
  vars_map["mindphin_metjet"] = var("mindphin_metjet","#Delta#phi_{N}", 10, 0., 40.);

  //---------------- DEFINE VARIABLE PAIRS 2D HISTOGRAMS ----------------
  vector<pair<TString,TString> > var_pairs; //use only variables defined in vars_map
  if (mode30){
    var_pairs.push_back(make_pair("ht30","mj_30"));
    var_pairs.push_back(make_pair("njets30","mj_30"));
    var_pairs.push_back(make_pair("nbl30","mj_30"));
    var_pairs.push_back(make_pair("met","mj_30"));
    var_pairs.push_back(make_pair("met","njets30"));
    var_pairs.push_back(make_pair("met","ht"));
    var_pairs.push_back(make_pair("mt","mj_30"));
    var_pairs.push_back(make_pair("mindphin_metjet","mj_30"));
  } else {
    var_pairs.push_back(make_pair("ht","mj_40"));
    var_pairs.push_back(make_pair("njets","mj_40"));
    var_pairs.push_back(make_pair("njets","ht"));
    var_pairs.push_back(make_pair("nbl","mj_40"));
    var_pairs.push_back(make_pair("met","mj_40"));
    var_pairs.push_back(make_pair("met","ht"));
    var_pairs.push_back(make_pair("met","njets"));
    var_pairs.push_back(make_pair("mt","mj_40"));
    var_pairs.push_back(make_pair("mindphin_metjet","mj_40"));
  }

  //---------------- DEFINE VARIABLE SETS 3D HISTOGRAMS ----------------
  vector<vector<TString> > var_sets(0,vector<TString>(0)); //use only variables defined in vars_map
  if (mode30) {
    var_sets.push_back(vector<TString>(0));
    var_sets[0].push_back("mj_30");   var_sets[0].push_back("njets30");   var_sets[0].push_back("met");
    var_sets.push_back(vector<TString>(0));
    var_sets[1].push_back("ht30");   var_sets[1].push_back("njets30");   var_sets[1].push_back("met");
  } else {
    var_sets.push_back(vector<TString>(0));
    var_sets[0].push_back("mj_40");   var_sets[0].push_back("njets");   var_sets[0].push_back("met");
    var_sets.push_back(vector<TString>(0));
    var_sets[1].push_back("ht");   var_sets[1].push_back("njets");   var_sets[1].push_back("met");
  }

  for (vector<sample>::iterator isamp = samples.begin(); isamp != samples.end(); isamp++) {

    TString outdir = "./out/jet40/mj_plots_"+isamp->name+".root";
    if (mode30) outdir= "./out/jet30/mj_plots_"+isamp->name+".root";
    TFile fout(outdir,"RECREATE");
    small_tree tree((basedir+isamp->filestr).Data());

    map<TString, h1d> h1d_map;
    map<TString, h2d> h2d_map;
    map<TString, h3d> h3d_map;

    //---------------- GENERATE 1D HISTOGRAMS FOR ALL VARIABLES ----------------
    for (map<TString, var>::iterator ivar = vars_map.begin(); ivar!=vars_map.end(); ivar++){
      for (vector<seln>::iterator iseln = selns.begin(); iseln!=selns.end(); iseln++){
        TString hname = iseln->name+"_"+ivar->second.name+"_"+isamp->name;
        h1d_map[hname] = h1d(ivar->second,hname);
      }
    }

    //---------------- GENERATE 2D HISTOGRAMS FOR ALL VARIABLE PAIRS ----------------
    for (vector<pair<TString,TString> >::iterator ipair=var_pairs.begin(); ipair!=var_pairs.end(); ipair++){
      for (vector<seln>::iterator iseln = selns.begin(); iseln!=selns.end(); iseln++){
        TString hname = iseln->name+"_"+ipair->first+"_"+ipair->second+"_"+isamp->name;
        h2d_map[hname] = h2d(vars_map[ipair->first], vars_map[ipair->second], hname);
      }
    }

    //---------------- GENERATE 3D HISTOGRAMS FOR ALL VARIABLE SETS ----------------
    for (vector<vector<TString> >::iterator iset=var_sets.begin(); iset!=var_sets.end(); iset++){
      for (vector<seln>::iterator iseln = selns.begin(); iseln!=selns.end(); iseln++){
        TString hname = iseln->name+"_"+iset->at(0)+"_"+iset->at(1)+"_"+iset->at(2)+"_"+isamp->name;
        h3d_map[hname] = h3d(vars_map[iset->at(0)], vars_map[iset->at(1)], vars_map[iset->at(2)], hname);
      }
    }

    //---------------- EVENT LOOP ----------------
    const long nent = tree.GetEntries();
    cout<<"Number of events to run over: "<<nent<<endl;
    for (long ientry=0; ientry<nent; ientry++){
      if (ientry%100000==0) cout<<"Processed events: "<<ientry<<endl;
      tree.GetEntry(ientry);
      double weight = tree.weight()*5;

      map<TString, double> flt_val;
      // for the sake of filling histos all are floats...
      flt_val["mj_30"] = tree.mj_30();
      flt_val["ht30"] = tree.ht30();
      flt_val["njets30"] = double(tree.njets30())+0.5; 
      flt_val["nbl30"] = double(tree.nbl30())+0.5; 
      flt_val["mj_30"] = tree.mj_30();
      flt_val["ht"] = tree.ht();
      flt_val["njets"] = double(tree.njets())+0.5; 
      flt_val["nbl"] = double(tree.nbl())+0.5; 
      flt_val["mj_40"] = tree.mj_40();
      flt_val["met"] = tree.met(); 
      flt_val["mt"] = tree.mt(); 
      flt_val["mindphin_metjet"] = tree.mindphin_metjet(); 
      flt_val["nfjets_40"] = tree.nfjets_40(); 

      double mj_40_fjm70 = 0.;
      for (vector<float>::iterator fjmass=tree.fjets_40_m().begin(); fjmass!=tree.fjets_40_m().end(); fjmass++){
        if (*fjmass>70.) mj_40_fjm70 += *fjmass;
      }

      flt_val["mj_40_fjm70"] = mj_40_fjm70;


      for (vector<seln>::iterator iseln = selns.begin(); iseln!=selns.end(); iseln++){
      
        if (!passSelection(tree, *iseln)) continue;

        for (map<TString, var>::iterator ivar = vars_map.begin(); ivar!=vars_map.end(); ivar++){
          TString hname = iseln->name+"_"+ivar->second.name+"_"+isamp->name;
          if (ivar->second.name!="fjets_40_m") {
            h1d_map[hname].fill(flt_val[ivar->second.name], weight);
          } else {
            for (vector<float>::iterator fjmass=tree.fjets_40_m().begin(); fjmass!=tree.fjets_40_m().end(); fjmass++)
              h1d_map[hname].fill(*fjmass, weight);
          }
        }

        for (vector<pair<TString,TString> >::iterator ipair=var_pairs.begin(); ipair!=var_pairs.end(); ipair++){
          TString hname = iseln->name+"_"+ipair->first+"_"+ipair->second+"_"+isamp->name;
          h2d_map[hname].fill(flt_val[ipair->first], flt_val[ipair->second], weight);
        }

        for (vector<vector<TString> >::iterator iset=var_sets.begin(); iset!=var_sets.end(); iset++){
          TString hname = iseln->name+"_"+iset->at(0)+"_"+iset->at(1)+"_"+iset->at(2)+"_"+isamp->name;
          h3d_map[hname].fill(flt_val[iset->at(0)], flt_val[iset->at(1)], flt_val[iset->at(2)], weight);
        }

      }    
    }

  fout.Write();
  fout.Close();
  }

  return 0;
}

bool passIsolation(const small_tree &tree, int ilep, bool isElectron, bool isveto){

  if (isElectron) {
    if (isveto)
      return ((tree.els_miniso_tr10()[ilep] < 0.1) || (tree.els_reliso_r02()[ilep] < 0.1));
    else 
      return ((tree.els_miniso_tr10()[ilep] < 0.1) || (tree.els_reliso_r02()[ilep] < 0.1));
  } else { 
    if (isveto)
      return ((tree.mus_miniso_tr10()[ilep] < 0.4) || (tree.mus_reliso_r02()[ilep] < 0.4));
    else 
      return ((tree.mus_miniso_tr10()[ilep] < 0.4) || (tree.mus_reliso_r02()[ilep] < 0.4));
  } 
}

bool passSelection(const small_tree &tree, const seln &iseln){

  size_t nveto_mus(0), nveto_els(0);
  vector<int> sigels_index = vector<int>(0);
  vector<int> sigmus_index = vector<int>(0);
  //---------- ELECTRONS ----------------      
  for (size_t iel=0; iel<tree.els_pt().size(); iel++){
    if (!tree.els_ispf()[iel]) continue;
    if (tree.els_pt()[iel] < iseln.veto_lep_pt_min) continue;
    if (passIsolation(tree, iel, /*isElectron*/ true, /*isveto*/ true)) nveto_els++;


    if (!(tree.els_sigid()[iel])) continue;
    if (tree.els_pt()[iel] < iseln.sig_lep_pt_min) continue;
    if (passIsolation(tree, iel, /*isElectron*/ true, /*isveto*/ false)) sigels_index.push_back(iel);
  }
  //---------- MUONS ----------------
  for (size_t imu=0; imu<tree.mus_pt().size(); imu++){
    if (tree.mus_pt()[imu] < iseln.veto_lep_pt_min) continue;
    if (passIsolation(tree, imu, /*isElectron*/ false, /*isveto*/ true)) nveto_mus++;

    if (!(tree.mus_sigid()[imu])) continue;
    if (tree.mus_pt()[imu] < iseln.sig_lep_pt_min) continue;
    if (passIsolation(tree, imu, /*isElectron*/ false, /*isveto*/ false)) sigmus_index.push_back(imu);
  }

  if ((sigels_index.size()+sigmus_index.size())!=iseln.nleps || (nveto_els+nveto_mus)!=iseln.nleps) return false;

  if (iseln.nleps==1) {
    double mt = 0;
    if (sigels_index.size()>0) mt = sqrt(2*tree.els_pt()[sigels_index[0]]*tree.met()*(1-cos(tree.met_phi()-tree.els_phi()[sigels_index[0]])));
    else if (sigmus_index.size()>0) mt = sqrt(2*tree.mus_pt()[sigmus_index[0]]*tree.met()*(1-cos(tree.met_phi()-tree.mus_phi()[sigmus_index[0]])));
    if (tree.mt() <= iseln.mt_min || tree.mt() > iseln.mt_max) return false;
  }

  if (mode30) { 
    if (tree.ht30() <= iseln.ht_min || tree.ht30() > iseln.ht_max) return false;
    if (tree.mj_30() <= iseln.mj_min || tree.mj_30() > iseln.mj_max) return false;
    if (tree.njets30() <= iseln.njets_min || tree.njets30() > iseln.njets_max) return false;
    if (tree.nbl30() <= iseln.nbl_min || tree.nbl30() > iseln.nbl_max) return false;
  } else {
    if (tree.ht() <= iseln.ht_min || tree.ht() > iseln.ht_max) return false;
    if (tree.mj_40() <= iseln.mj_min || tree.mj_40() > iseln.mj_max) return false;
    if (tree.njets() <= iseln.njets_min || tree.njets() > iseln.njets_max) return false;
    if (tree.nbl() <= iseln.nbl_min || tree.nbl() > iseln.nbl_max) return false;
  }

  if (tree.met() <= iseln.met_min || tree.met() > iseln.met_max) return false;

  if (iseln.nleps==0 && tree.mindphin_metjet()<=4.0) return false;

  return true;
}

