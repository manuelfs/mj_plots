#include "mj_correl.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <set>
#include <map>
#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLorentzVector.h"

#include "RooStats/NumberCountingUtils.h"

#include "utilities.hpp"

using namespace std;

namespace{
  const double dblmax = numeric_limits<double>::max();
  const int intmax = numeric_limits<int>::max();
}

namespace msg { enum { dbg=0, info=1, warn=2, err=3}; }

int main(){

  unsigned msglvl = msg::info;

  bool mode30 = false;

  TString basedir = "/afs/cern.ch/user/m/manuelf/work/ucsb/15-01-30_skim/";
  if (mode30) basedir = "/afs/cern.ch/user/m/manuelf/work/ucsb/15-01-30_skim/ht30/";

  vector<sample> samples;
  samples.push_back(sample("*_TTJets*", "ttbar", true));
  samples.push_back(sample("*-T1tttt_2J_mGl-1500_mLSP-100_*PU20*", "T1tttt1500", true));
  samples.push_back(sample("*-T1tttt_2J_mGl-1200_mLSP-800_*PU20*", "T1tttt1200", true));
  samples.push_back(sample("*-T1qqqq_2J_mGl-1000_mLSP-800_*", "T1qqqq1000", false));
  samples.push_back(sample("*-T1qqqq_2J_mGl-1400_mLSP-100_*", "T1qqqq1400", false));
  // samples.push_back(sample("*-T2tt_2J_mStop-650_mLSP-325_*", "T2tt650", true)); skim not available
  // samples.push_back(sample("*-T2tt_2J_mStop-850_mLSP-100_*", "T2tt850", true));
  samples.push_back(sample("*_WJetsToLNu*", "wjets", false));
  samples.push_back(sample("*_QCD_*", "qcd", false));

  //---------------- DEFINE EVENT SELECTIONS ----------------
  // double els_iso_max = 0.1;
  // double mus_iso_max = 0.4;
  
  vector<seln> selns;
//  selns.push_back(seln("0L-base",
//                /*nleps_*/ 0, /*nisotks_max_*/ intmax,
//                /*sig_lep_pt_min*/ 15.,/*veto_lep_pt_min*/ 15.,
//                /*ht_min_*/ 500., /*ht_max_*/ dblmax,
//                /*mj_min_*/ 0., /*mj_max_*/ dblmax,
//                /*met_min_*/ 200., /*met_max_*/ dblmax,
//                /*mt_min_*/ 0, /*mt_max_*/ dblmax,
//                /*njets_min_*/ 0, /*njets_max_*/ intmax,
//                /*nbl_min_*/ 0, /*nbl_max_*/ intmax));
//  selns.push_back(seln("1L-base",
//                  /*nleps_*/ 1, /*nisotks_max_*/ intmax,
//                  /*sig_lep_pt_min*/ 15.,/*veto_lep_pt_min*/ 15.,
//                  /*ht_min_*/ 500., /*ht_max_*/ dblmax,
//                  /*mj_min_*/ 0., /*mj_max_*/ dblmax,
//                  /*met_min_*/ 200., /*met_max_*/ dblmax,
//                  /*mt_min_*/ 0, /*mt_max_*/ dblmax,
//                  /*njets_min_*/ 0, /*njets_max_*/ intmax,
//                  /*nbl_min_*/ 0, /*nbl_max_*/ intmax));
//
//  selns.push_back(seln("2L-base",
//                /*nleps_*/ 2, /*nisotks_max_*/ intmax,
//                /*sig_lep_pt_min*/ 15.,/*veto_lep_pt_min*/ 15.,
//                /*ht_min_*/ 500., /*ht_max_*/ dblmax,
//                /*mj_min_*/ 0., /*mj_max_*/ dblmax,
//                /*met_min_*/ 200., /*met_max_*/ dblmax,
//                /*mt_min_*/ 0, /*mt_max_*/ dblmax,
//                /*njets_min_*/ 0, /*njets_max_*/ intmax,
//                /*nbl_min_*/ 0, /*nbl_max_*/ intmax));

    selns.push_back(seln("0L-sign",
                /*nleps_*/ 0, /*nisotks_max_*/ intmax,
                /*sig_lep_pt_min*/ 15.,/*veto_lep_pt_min*/ 15.,
                /*ht_min_*/ 500., /*ht_max_*/ dblmax,
                /*mj_min_*/ 0., /*mj_max_*/ dblmax,
                /*met_min_*/ 200., /*met_max_*/ dblmax,
                /*mt_min_*/ 0, /*mt_max_*/ dblmax,
                /*njets_min_*/ 6, /*njets_max_*/ intmax,
                /*nbl_min_*/ 2, /*nbl_max_*/ intmax));
  selns.push_back(seln("1L-sign",
                  /*nleps_*/ 1, /*nisotks_max_*/ intmax,
                  /*sig_lep_pt_min*/ 15.,/*veto_lep_pt_min*/ 15.,
                  /*ht_min_*/ 500., /*ht_max_*/ dblmax,
                  /*mj_min_*/ 0., /*mj_max_*/ dblmax,
                  /*met_min_*/ 200., /*met_max_*/ dblmax,
                  /*mt_min_*/ 0, /*mt_max_*/ dblmax,
                  /*njets_min_*/ 6, /*njets_max_*/ intmax,
                  /*nbl_min_*/ 2, /*nbl_max_*/ intmax));

  selns.push_back(seln("2L-sign",
                /*nleps_*/ 2, /*nisotks_max_*/ intmax,
                /*sig_lep_pt_min*/ 15.,/*veto_lep_pt_min*/ 15.,
                /*ht_min_*/ 500., /*ht_max_*/ dblmax,
                /*mj_min_*/ 0., /*mj_max_*/ dblmax,
                /*met_min_*/ 200., /*met_max_*/ dblmax,
                /*mt_min_*/ 0, /*mt_max_*/ dblmax,
                /*njets_min_*/ 6, /*njets_max_*/ intmax,
                /*nbl_min_*/ 2, /*nbl_max_*/ intmax));

  msgsvc(msglvl, msg::dbg,"Defined selections...");

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
    vars_map["lead_fjets_40_m"] = var("lead_fjets_40_m", "lead-m large-R jet mass [GeV]", 30, 0., 1500.);
    vars_map["lead_fjets_40_pt"] = var("lead_fjets_40_pt", "lead-m large-R jet p_{T} [GeV]", 30, 0., 1500.);
    vars_map["nfjets_40"] = var("nfjets_40", "large-R jet multiplicity", 20, 0., 20.);
  }
  vars_map["met"] = var("met", "MET [GeV]", 30, 0., 1500.);
  vars_map["mt"] = var("mt", "m_{T} [GeV]", 20, 0., 1000.);
  vars_map["mindphin_metjet"] = var("mindphin_metjet","#Delta#phi_{N}", 10, 0., 40.);
  vars_map["pt_tt"] = var("pt_tt", "p^{gen}_{T}(t#bar{t}) [GeV]", 75, 0., 1500.);
  vars_map["lead_pt_top"] = var("lead_pt_top", "leading top p^{gen}_{T} [GeV]", 30, 0., 1500.);
  vars_map["sublead_pt_top"] = var("sublead_pt_top", "sub-leading top p^{gen}_{T} [GeV]", 30, 0., 1500.);
  vars_map["lead_pt_gluon"] = var("lead_pt_gluon", "leading gluon p^{gen}_{T} [GeV]", 30, 0., 1500.);

  msgsvc(msglvl, msg::dbg,"Defined variables...");

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
  } else {
    var_pairs.push_back(make_pair("ht","mj_40"));
    var_pairs.push_back(make_pair("njets","mj_40"));
    var_pairs.push_back(make_pair("njets","ht"));
    var_pairs.push_back(make_pair("nbl","mj_40"));
    var_pairs.push_back(make_pair("met","mj_40"));
    var_pairs.push_back(make_pair("met","ht"));
    var_pairs.push_back(make_pair("met","njets"));
    var_pairs.push_back(make_pair("pt_tt","lead_pt_top"));
    var_pairs.push_back(make_pair("pt_tt","lead_pt_gluon"));
    var_pairs.push_back(make_pair("pt_tt","mj_40"));
    var_pairs.push_back(make_pair("lead_pt_top","mj_40"));
    var_pairs.push_back(make_pair("pt_tt","ht"));
    var_pairs.push_back(make_pair("lead_pt_top","ht"));
  }

  msgsvc(msglvl, msg::dbg,"Defined correlations...");
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

  msgsvc(msglvl, msg::dbg,"Defined 3D plots...");

  for (vector<sample>::iterator isamp = samples.begin(); isamp != samples.end(); isamp++) {

    msgsvc(msglvl, msg::dbg, "Running over sample: " + isamp->name);
    TString outdir = "./out/sign/mj_plots_"+isamp->name+".root";
    TFile fout(outdir,"RECREATE");
    small_tree tree((basedir+isamp->filestr).Data());

    map<TString, h1d> h1d_map;
    map<TString, h2d> h2d_map;
    map<TString, h3d> h3d_map;

    //---------------- CREATE 1D HISTOGRAMS FOR ALL VARIABLES ----------------
    for (map<TString, var>::iterator ivar = vars_map.begin(); ivar!=vars_map.end(); ivar++){
      for (vector<seln>::iterator iseln = selns.begin(); iseln!=selns.end(); iseln++){
        TString hname = iseln->name+"_"+ivar->second.name+"_"+isamp->name;
        msgsvc(msglvl, msg::dbg, "Creating histogram: " + hname);
        h1d_map[hname] = h1d(ivar->second,hname);
      }
    }

    //---------------- CREATE 2D HISTOGRAMS FOR ALL VARIABLE PAIRS ----------------
    for (vector<pair<TString,TString> >::iterator ipair=var_pairs.begin(); ipair!=var_pairs.end(); ipair++){
      for (vector<seln>::iterator iseln = selns.begin(); iseln!=selns.end(); iseln++){
        TString hname = iseln->name+"_"+ipair->first+"_"+ipair->second+"_"+isamp->name;
        msgsvc(msglvl, msg::dbg, "Creating histogram: " + hname);
        h2d_map[hname] = h2d(vars_map[ipair->first], vars_map[ipair->second], hname);
      }
    }

    //---------------- CREATE 3D HISTOGRAMS FOR ALL VARIABLE SETS ----------------
    for (vector<vector<TString> >::iterator iset=var_sets.begin(); iset!=var_sets.end(); iset++){
      for (vector<seln>::iterator iseln = selns.begin(); iseln!=selns.end(); iseln++){
        TString hname = iseln->name+"_"+iset->at(0)+"_"+iset->at(1)+"_"+iset->at(2)+"_"+isamp->name;
        msgsvc(msglvl, msg::dbg, "Creating histogram: " + hname);
        h3d_map[hname] = h3d(vars_map[iset->at(0)], vars_map[iset->at(1)], vars_map[iset->at(2)], hname);
      }
    }

    //---------------- EVENT LOOP ----------------
    const size_t nent = tree.GetEntries();
    msgsvc(msglvl, msg::info, TString::Format("Number of events to run over: %u", unsigned(nent)));
    for (size_t ientry(0); ientry<nent; ientry++){
      if (ientry%100000==0) cout<<"INFO:: Processed events: "<<ientry<<endl;
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

      // ------------ CALCULATE CUSTOM VARIABLES --------------
      
      // ------ MJ with massive fat jets and leading mass jet vars ----------
      flt_val["mj_40_fjm70"] = 0.;
      flt_val["lead_fjets_40_m"] = 0.; 
      flt_val["lead_fjets_40_pt"] = 0.; 
      for (size_t i(0); i<tree.fjets_40_m().size(); i++){
        if (tree.fjets_40_m()[i]>70.) flt_val["mj_40_fjm70"] += tree.fjets_40_m()[i];
        if (tree.fjets_40_m()[i]>flt_val["lead_fjets_40_m"]){
          flt_val["lead_fjets_40_m"] = tree.fjets_40_m()[i]; 
          flt_val["lead_fjets_40_pt"] = tree.fjets_40_pt()[i]; 
        }
      }

      // ------ Generator level top pt ----------
      flt_val["lead_pt_top"] = -1.*dblmax;
      flt_val["sublead_pt_top"] = -1.*dblmax;
      flt_val["lead_pt_gluon"] = -1.*dblmax;
      vector<int> tops; 
      vector<double> tops_pt; 
      for (size_t i(0); i<tree.mc_id().size(); i++){
        if (abs(tree.mc_id()[i])==6) {
          tops.push_back(i);
          tops_pt.push_back(tree.mc_pt()[i]);
        } else if (abs(tree.mc_id()[i])==21) {
          if (tree.mc_pt()[i] > flt_val["lead_pt_gluon"]) flt_val["lead_pt_gluon"] = tree.mc_pt()[i];
        }
      }
      sort(tops_pt.begin(), tops_pt.end());
      size_t ntops = tops_pt.size();
      if (ntops>0) flt_val["lead_pt_top"] = tops_pt[ntops-1];
      if (ntops>1) flt_val["sublead_pt_top"] = tops_pt[ntops-2];

      // ------ Generator level pt of the ttbar system ----------
      if (isamp->name=="ttbar"){
        if (ntops!=2){
          cout<<"WARNING:: Bad truth record, number of tops = "<<ntops<<endl;
          flt_val["pt_tt"] = dblmax; 
        } else {
          flt_val["pt_tt"] = sqrt(pow(tree.mc_pt()[tops[0]],2) + pow(tree.mc_pt()[tops[1]],2) +
                  2*tree.mc_pt()[tops[0]]*tree.mc_pt()[tops[1]]*cos(DeltaPhi(tree.mc_phi()[tops[0]],tree.mc_phi()[tops[1]])));
        }
      } else {
        flt_val["pt_tt"] = dblmax;
      }


      for (vector<seln>::iterator iseln = selns.begin(); iseln!=selns.end(); iseln++){
      
        if (!passSelection(tree, *iseln, mode30)) continue;

        for (map<TString, var>::iterator ivar = vars_map.begin(); ivar!=vars_map.end(); ivar++){
          TString hname = iseln->name+"_"+ivar->second.name+"_"+isamp->name;
          h1d_map[hname].fill(flt_val[ivar->second.name], weight);
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

bool passSelection(const small_tree &tree, const seln &iseln, bool mode30){

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
    if (tree.njets30() < iseln.njets_min || tree.njets30() > iseln.njets_max) return false;
    if (tree.nbl30() < iseln.nbl_min || tree.nbl30() > iseln.nbl_max) return false;
  } else {
    if (tree.ht() <= iseln.ht_min || tree.ht() > iseln.ht_max) return false;
    if (tree.mj_40() <= iseln.mj_min || tree.mj_40() > iseln.mj_max) return false;
    if (tree.njets() < iseln.njets_min || tree.njets() > iseln.njets_max) return false;
    if (tree.nbl() < iseln.nbl_min || tree.nbl() > iseln.nbl_max) return false;
  }

  if (tree.met() <= iseln.met_min || tree.met() > iseln.met_max) return false;

  if (iseln.nleps==0 && tree.mindphin_metjet()<=4.0) return false;

  return true;
}

void msgsvc(const unsigned &setlvl, const unsigned &lvl, const TString &mymsg) {
  if (lvl>=setlvl) {
    if(lvl == msg::dbg) cout<<"DEBUG:: "<<mymsg<<endl;
    else if(lvl == msg::info) cout<<"INFO:: "<<mymsg<<endl;
    else if(lvl == msg::warn) cout<<"WARNING:: "<<mymsg<<endl;
    else if(lvl == msg::err) cout<<"ERROR:: "<<mymsg<<endl;

  }
}
