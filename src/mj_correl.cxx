#include "mj_correl.hpp"
#include <iostream>
#include <iomanip>
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
#include "TVector2.h"
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

  unsigned msglvl = msg::dbg;
  // unsigned msglvl = msg::info;

  // TString basedir = "/cms5r0/manuelf/root/archive/15-01-30/skim/";
  TString basedir = "/cms5r0/ald77/archive/current/";

  vector<sample> samples; 
  samples.push_back(sample("*_TTJets*", "ttbar", true));
  // samples.push_back(sample("*_TTZ*", "ttZ", true));
  // samples.push_back(sample("*-T1tttt_2J_mGl-1500_mLSP-100_*PU20*", "T1tttt1500", true));
  // samples.push_back(sample("*-T1tttt_2J_mGl-1200_mLSP-800_*PU20*", "T1tttt1200", true));
  // samples.push_back(sample("*-T1qqqq_2J_mGl-1000_mLSP-800_*", "T1qqqq1000", false));
  // samples.push_back(sample("*-T1qqqq_2J_mGl-1400_mLSP-100_*", "T1qqqq1400", false));
  // samples.push_back(sample("*-T2tt_2J_mStop-650_mLSP-325_*", "T2tt650", true)); skim not available
  // samples.push_back(sample("*-T2tt_2J_mStop-850_mLSP-100_*", "T2tt850", true));
  // samples.push_back(sample("*_WJetsToLNu*", "wjets", false));
  // samples.push_back(sample("*_QCD_*", "qcd", false));

  //---------------- DEFINE EVENT SELECTIONS ----------------
  vector<seln> selns;
  selns.push_back(seln("nl1-ht750-met250-mt150-nj6-nb2", /*nleps_*/ 1,  intmax, /*sig_pt*/ 20.,/*veto_pt*/ 10., /*ht*/ 750., dblmax, /*mj*/ 0.,  dblmax, /*met*/ 250.,  dblmax, /*mt*/ 150,  dblmax, /*njets*/ 6,  intmax, /*nbl*/ 2,  intmax));
  selns.push_back(seln("nl1-ht500-met200-mt150-nj0-nb0", /*nleps_*/ 1,  intmax, /*sig_pt*/ 20.,/*veto_pt*/ 10., /*ht*/ 500., dblmax, /*mj*/ 0.,  dblmax, /*met*/ 200.,  dblmax, /*mt*/ 150,  dblmax, /*njets*/0,  intmax, /*nbl*/ 0,  intmax));
  selns.push_back(seln("nl1-ht500-met200-mt150-nj6-nb0", /*nleps_*/ 1,  intmax, /*sig_pt*/ 20.,/*veto_pt*/ 10., /*ht*/ 500., dblmax, /*mj*/ 0.,  dblmax, /*met*/ 200.,  dblmax, /*mt*/ 150,  dblmax, /*njets*/6,  intmax, /*nbl*/ 0,  intmax));
  selns.push_back(seln("nl1-ht500-met200-mt0-nj6-nb0", /*nleps_*/ 1,  intmax, /*sig_pt*/ 20.,/*veto_pt*/ 10., /*ht*/ 500., dblmax, /*mj*/ 0.,  dblmax, /*met*/ 200.,  dblmax, /*mt*/ 0,  dblmax, /*njets*/6,  intmax, /*nbl*/ 0,  intmax));
  selns.push_back(seln("nl1-ht500-met200-mt0-nj6-nb0-mj0to200", /*nleps_*/ 1,  intmax, /*sig_pt*/ 20.,/*veto_pt*/ 10., /*ht*/ 500., dblmax, /*mj*/ 0.,  200., /*met*/ 200.,  dblmax, /*mt*/ 0,  dblmax, /*njets*/6,  intmax, /*nbl*/ 0,  intmax));
  selns.push_back(seln("nl1-ht500-met200-mt0-nj6-nb0-mj200to400", /*nleps_*/ 1,  intmax, /*sig_pt*/ 20.,/*veto_pt*/ 10., /*ht*/ 500., dblmax, /*mj*/ 200.,  400., /*met*/ 200.,  dblmax, /*mt*/ 0,  dblmax, /*njets*/6,  intmax, /*nbl*/ 0,  intmax));
  selns.push_back(seln("nl1-ht500-met200-mt0-nj6-nb0-mj400to600", /*nleps_*/ 1,  intmax, /*sig_pt*/ 20.,/*veto_pt*/ 10., /*ht*/ 500., dblmax, /*mj*/ 400.,  600., /*met*/ 200.,  dblmax, /*mt*/ 0,  dblmax, /*njets*/6,  intmax, /*nbl*/ 0,  intmax));
  selns.push_back(seln("nl1-ht500-met200-mt0-nj6-nb0-mj600toInf", /*nleps_*/ 1,  intmax, /*sig_pt*/ 20.,/*veto_pt*/ 10., /*ht*/ 500., dblmax, /*mj*/ 600.,  dblmax, /*met*/ 200.,  dblmax, /*mt*/ 0,  dblmax, /*njets*/6,  intmax, /*nbl*/ 0,  intmax));
  msgsvc(msglvl, msg::dbg,"Defined selections...");

  //---------------- DEFINE VARIABLES FOR ALL PLOTS ----------------
  map<TString, var> vars_map;
  vars_map["mj"] = var("mj", "Sum(m_{J}) [GeV]", 30, 0., 1500.);
  vars_map["ht"] = var("ht", "H_{T} [GeV]", 35, 500., 4000.);
  vars_map["njets"] = var("njets", "jet multiplicity", 20, 0., 20.);
  vars_map["nbl"] = var("nbl", "b-jet multiplicity", 6, 0., 6.);
  vars_map["fjm1"] = var("fjm1", "m_{j1} [GeV]", 10, 0., 500.);
  vars_map["fjm2"] = var("fjm2", "m_{j2} [GeV]", 10, 0., 500.);
  vars_map["leadfjpt"] = var("leadfjpt", "lead-m large-R jet p_{T} [GeV]", 30, 0., 1500.);
  vars_map["nfjets"] = var("nfjets", "large-R jet multiplicity", 20, 0., 20.);
  vars_map["met"] = var("met", "MET [GeV]", 30, 0., 1500.);
  vars_map["mt"] = var("mt", "m_{T} [GeV]", 20, 0., 1000.);
  vars_map["isr"] = var("isr", "ISR [GeV]", 10, 0., 1000.);
  vars_map["isr1pt"] = var("isr1pt", "ISR 1st parton [GeV]", 10, 0., 1000.);
  vars_map["isr2pt"] = var("isr2pt", "ISR 2nd parton [GeV]", 10, 0., 1000.);
  vars_map["isr3pt"] = var("isr3pt", "ISR 3rd parton [GeV]", 10, 0., 1000.);
  vars_map["nisrjets"] = var("nisrjets", "ISR parton multiplicity", 6, 0., 6.);
  vars_map["ptt"] = var("ptt", "p^{gen}_{T}(t#bar{t}) [GeV]", 10, 0., 1000.);
  vars_map["avetoppt"] = var("avetoppt", "(p^{gen}_{T}(t)+p^{gen}_{T}(#bar{t}))/2 [GeV]", 10, 0., 1000.);
  vars_map["dphi_tt"] = var("dphi_tt", "#Delta#phi^{gen}(t,#bar{t})", 10, 0., 3.14);
  vars_map["dphi_fjm1_fjm2"] = var("dphi_fjm1_fjm2", "#Delta#phi^{gen}(m_{j1},m_{j2})", 10, 0., 3.14);
  vars_map["mindphi_isr_top"] = var("mindphi_isr_top", "min(#Delta#phi(t_{i},ISR_j))", 10, 0., 3.14);
  vars_map["lead_pt_top"] = var("lead_pt_top", "leading top p^{gen}_{T} [GeV]", 30, 0., 1500.);
  vars_map["sublead_pt_top"] = var("sublead_pt_top", "sub-leading top p^{gen}_{T} [GeV]", 30, 0., 1500.);
  vars_map["leadptglu"] = var("leadptglu", "leading gluon p^{gen}_{T} [GeV]", 10, 0., 1000.);
  msgsvc(msglvl, msg::dbg,"Defined variables...");

  //---------------- DEFINE VARIABLE PAIRS 2D HISTOGRAMS ----------------
  vector<pair<TString,TString> > var_pairs; //use only variables defined in vars_map
  var_pairs.push_back(make_pair("ht","mj"));
  // var_pairs.push_back(make_pair("njets","mj"));
  // var_pairs.push_back(make_pair("njets","ht"));
  // var_pairs.push_back(make_pair("nbl","mj"));
  // var_pairs.push_back(make_pair("met","mj"));
  // var_pairs.push_back(make_pair("met","ht"));
  // var_pairs.push_back(make_pair("met","njets"));
  var_pairs.push_back(make_pair("lead_pt_top", "ptt"));
  var_pairs.push_back(make_pair("lead_pt_top", "leadfjpt"));
  var_pairs.push_back(make_pair("isr1pt", "avetoppt"));
  var_pairs.push_back(make_pair("isr2pt", "avetoppt"));
  var_pairs.push_back(make_pair("isr3pt", "avetoppt"));
  var_pairs.push_back(make_pair("ptt","mj"));
  var_pairs.push_back(make_pair("lead_pt_top","mj"));
  var_pairs.push_back(make_pair("dphi_tt","ptt"));
  var_pairs.push_back(make_pair("dphi_tt","avetoppt"));
  var_pairs.push_back(make_pair("dphi_tt","mj"));
  var_pairs.push_back(make_pair("dphi_tt","dphi_fjm1_fjm2"));
  var_pairs.push_back(make_pair("nisrjets","mj"));
  var_pairs.push_back(make_pair("nisrjets","ptt"));
  var_pairs.push_back(make_pair("mindphi_isr_top","mj"));
  var_pairs.push_back(make_pair("mindphi_isr_top","fjm1"));
  var_pairs.push_back(make_pair("dphi_tt","nisrjets"));
  var_pairs.push_back(make_pair("lead_pt_top","dphi_tt"));
  var_pairs.push_back(make_pair("fjm1","dphi_tt"));
  msgsvc(msglvl, msg::dbg,"Defined correlations...");

  //---------------- DEFINE VARIABLE SETS 3D HISTOGRAMS ----------------
  vector<vector<TString> > var_sets(0,vector<TString>(0)); //use only variables defined in vars_map
  int ivarset = 0;
  var_sets.push_back(vector<TString>(0));
  var_sets[ivarset].push_back("ptt");   var_sets[ivarset].push_back("dphi_tt");   var_sets[ivarset].push_back("mj");

  // ivarset++;
  // var_sets.push_back(vector<TString>(0));
  // var_sets[ivarset].push_back("ht");   var_sets[ivarset].push_back("njets");   var_sets[ivarset].push_back("met");

  // ivarset++;
  // var_sets.push_back(vector<TString>(0));
  // var_sets[ivarset].push_back("mj");   var_sets[ivarset].push_back("njets");   var_sets[ivarset].push_back("met");

  msgsvc(msglvl, msg::dbg,"Defined 3D plots...");

  //---------------- LOOP OVER SAMPLES ----------------
  for (vector<sample>::iterator isamp = samples.begin(); isamp != samples.end(); isamp++) {

    msgsvc(msglvl, msg::dbg, "Running over sample: " + isamp->name);
    TString outdir = "./ntup/mj_plots_"+isamp->name+".root";
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

    //---------------- LOOP OVER EVENTS ----------------
    const size_t nent = tree.GetEntries();
    msgsvc(msglvl, msg::info, TString::Format("Number of events to run over: %u", unsigned(nent)));
    for (size_t ientry(0); ientry<nent; ientry++){
      if (ientry%100000==0) cout<<"INFO:: Processed events: "<<ientry<<endl;

      // if ((tree.mc_type()&0x0F00)>=0x200) continue;
      // cout<<"Event passed"<<endl;
      tree.GetEntry(ientry);
      double weight = tree.weight()*5;

      map<TString, double> flt_val;
      // for the sake of filling histos all are floats...
      flt_val["ht"] = tree.ht();
      flt_val["njets"] = double(tree.njets())+0.5; 
      flt_val["nbl"] = double(tree.nbl())+0.5; 
      flt_val["mj"] = tree.mj();
      flt_val["met"] = tree.met(); 
      flt_val["mt"] = tree.mt(); 
      flt_val["nfjets"] = tree.nfjets(); 

      // ------------ CALCULATE CUSTOM VARIABLES --------------
      
      // ------ MJ with massive fat jets and leading mass jet vars ----------
      flt_val["fjm1"] = 0.; 
      flt_val["leadfjpt"] = 0.; 
      vector<TLorentzVector> fjs;
      for (size_t i(0); i<tree.fjets_m().size(); i++){
        TLorentzVector ifj;
        ifj.SetPtEtaPhiM(tree.fjets_pt()[i],tree.fjets_eta()[i],tree.fjets_phi()[i],tree.fjets_m()[i]);; 
        fjs.push_back(ifj);        
      }
      sort(fjs.begin(), fjs.end(), comp_mass);
      flt_val["dphi_fjm1_fjm2"] = fjs[0].DeltaPhi(fjs[1]);
      
      msgsvc(msglvl, msg::dbg, "Looking at the generator level ");
      // ------ Generator level ----------
      flt_val["lead_pt_top"] = -1.*dblmax;
      flt_val["sublead_pt_top"] = -1.*dblmax;
      flt_val["leadptglu"] = -1.*dblmax;
      flt_val["nisrjets"] = 0;
      TVector2 isr(1e-5,1e-5);
      vector<int> tops; 
      vector<double> tops_pt; 
      vector<double> isrs_pt; 
      vector<double> tops_phi; 
      vector<double> isrs_phi; 
      // cout<<"==========================================="<<endl;
      for (size_t i(0); i<tree.mc_id().size(); i++){
        // cout<<setw(10)<<tree.mc_id()[i]<<setw(10)<<tree.mc_mom()[i]<<endl;
            // cout<<"Id: "<<tree.mc_id()[i]<<"   Mom "<<tree.mc_mom()[i]<<endl;
        if (abs(tree.mc_id()[i])==6 && tree.mc_status()[i]==22) {
          tops.push_back(i);
          tops_pt.push_back(tree.mc_pt()[i]);
          tops_phi.push_back(tree.mc_phi()[i]);
        } else {
          if (abs(tree.mc_id()[i])==21 && tree.mc_pt()[i] > flt_val["leadptglu"]) flt_val["leadptglu"] = tree.mc_pt()[i];
          if (tree.mc_mom()[i]>100) {
            isrs_pt.push_back(tree.mc_pt()[i]);
            isrs_phi.push_back(tree.mc_phi()[i]);
            TVector2 ptcl; ptcl.SetMagPhi(tree.mc_pt()[i], tree.mc_phi()[i]);
            isr = isr + ptcl;
            flt_val["nisrjets"] = flt_val["nisrjets"]+1.;
          }
        }
      }
      flt_val["isr"] = isr.Mod();      

      sort(isrs_pt.begin(), isrs_pt.end());
      if (isrs_pt.size()>0) flt_val["isr1pt"] = isrs_pt[0];
      if (isrs_pt.size()>1) flt_val["isr2pt"] = isrs_pt[1];
      if (isrs_pt.size()>2) flt_val["isr3pt"] = isrs_pt[2];

      msgsvc(msglvl, msg::dbg, "Record top pTs ");
      sort(tops_pt.begin(), tops_pt.end());
      size_t ntops = tops_pt.size();
      if (ntops>0) flt_val["lead_pt_top"] = tops_pt[ntops-1];
      if (ntops>1) flt_val["sublead_pt_top"] = tops_pt[ntops-2];

      // min delta phi
      msgsvc(msglvl, msg::dbg, "Min delta phi between tops and isr ");
      double mindphi_isr_top = dblmax;
      for (size_t i(0); i<tops_phi.size(); i++){
        for (size_t j(0); j<isrs_phi.size(); j++){
          double dphi = DeltaPhi(tops_phi[i], isrs_phi[j]);
          if (dphi < mindphi_isr_top) mindphi_isr_top = dphi;
        }
      }
      flt_val["mindphi_isr_top"] = mindphi_isr_top;

      // ------ Generator level pt of the ttbar system ----------
      if (isamp->name=="ttbar"){
        if (ntops!=2){
          cout<<"WARNING:: Bad truth record, number of tops = "<<ntops<<endl;
          flt_val["ptt"] = dblmax; 
          flt_val["dphi_tt"] = dblmax;
          flt_val["avetoppt"] = dblmax;
        } else {
          TVector2 t1; t1.SetMagPhi(tree.mc_pt()[tops[0]], tree.mc_phi()[tops[0]]);
          TVector2 t2; t2.SetMagPhi(tree.mc_pt()[tops[1]], tree.mc_phi()[tops[1]]);
          TVector2 tt = t1+t2;
          flt_val["ptt"] = tt.Mod();
          flt_val["dphi_tt"] = t1.DeltaPhi(t2);
          flt_val["avetoppt"] = (tree.mc_pt()[tops[0]]+tree.mc_pt()[tops[1]])/2.;
        }
      } else {
        flt_val["ptt"] = dblmax;
        flt_val["dphi_tt"] = dblmax;
        flt_val["avetoppt"] = dblmax;
      }

      for (vector<seln>::iterator iseln = selns.begin(); iseln!=selns.end(); iseln++){
      
        if (!passSelection(tree, *iseln)) continue;

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

bool comp_mass(const TLorentzVector &fj1, const TLorentzVector &fj2){
  return (fj1.M()>fj2.M());
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

  if (tree.ht() <= iseln.ht_min || tree.ht() > iseln.ht_max) return false;
  if (tree.mj() <= iseln.mj_min || tree.mj() > iseln.mj_max) return false;
  if (tree.njets() < iseln.njets_min || tree.njets() > iseln.njets_max) return false;
  if (tree.nbm() < iseln.nbl_min || tree.nbm() > iseln.nbl_max) return false;
  

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
