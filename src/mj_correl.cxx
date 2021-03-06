#include "mj_correl.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <set>
#include <map>
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

namespace msg { enum {dbg=1, info=2, warn=3, err=4, truth=5}; }
// const unsigned msglvl = msg::truth;
unsigned msglvl = msg::info;
bool split_mctype = true;

int main(){

  // TString basedir = "/cms5r0/manuelf/root/archive/15-01-30/skim/";
  // TString basedir = "/cms5r0/ald77/archive/20150407/skims/";
  TString basedir = "/cms5r0/ald77/archive/20150414/";
  // TString basedir = "/cms5r0/ald77/archive/20150414/small_quick_TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola_Phys14DR-PU20bx25_PHYS14_25_V1-v1_MINIAODSIM_UCSB2290_v77_files50_batch11.root";
  // TString basedir = "/cms2r0/aovcharova/code/susy_cfa/out/";

  vector<TString> mctypes;
  // mctypes.push_back("ll");
  // mctypes.push_back("lt");
  // mctypes.push_back("lh");
  // mctypes.push_back("hh");
  mctypes.push_back("nl0");
  mctypes.push_back("nl1");
  mctypes.push_back("nl2");

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
  // selns.push_back(seln("nl1-ht500-met200-mt150-nj0-nb0", /*nleps_*/ 1,  intmax, /*sig_pt*/ 20.,/*veto_pt*/ 10., /*ht*/ 500., dblmax, /*mj*/ 0.,  dblmax, /*met*/ 200.,  dblmax, /*mt*/ 150,  dblmax, /*njets*/0,  intmax, /*nbl*/ 0,  intmax));
  // selns.push_back(seln("nl1-ht500-met200-mt150-nj6-nb0", /*nleps_*/ 1,  intmax, /*sig_pt*/ 20.,/*veto_pt*/ 10., /*ht*/ 500., dblmax, /*mj*/ 0.,  dblmax, /*met*/ 200.,  dblmax, /*mt*/ 150,  dblmax, /*njets*/6,  intmax, /*nbl*/ 0,  intmax));
  // selns.push_back(seln("nl1-ht500-met200-mt0-nj6-nb0", /*nleps_*/ 1,  intmax, /*sig_pt*/ 20.,/*veto_pt*/ 10., /*ht*/ 500., dblmax, /*mj*/ 0.,  dblmax, /*met*/ 200.,  dblmax, /*mt*/ 0,  dblmax, /*njets*/6,  intmax, /*nbl*/ 0,  intmax));
 
  selns.push_back(seln("nlna-ht500-met200-mt0-nj6-nb0", /*nleps_*/ -1,  intmax, /*sig_pt*/ 20.,/*veto_pt*/ 10., /*ht*/ 500., dblmax, /*mj*/ 0.,  dblmax, /*met*/ 200.,  dblmax, /*mt*/ 0,  dblmax, /*njets*/6,  intmax, /*nbl*/ 0,  intmax));
  // selns.push_back(seln("nl2-ht500-met200-mt0-nj6-nb0", /*nleps_*/ 2,  intmax, /*sig_pt*/ 20.,/*veto_pt*/ 10., /*ht*/ 500., dblmax, /*mj*/ 0.,  dblmax, /*met*/ 200.,  dblmax, /*mt*/ 0,  dblmax, /*njets*/6,  intmax, /*nbl*/ 0,  intmax));
  
  // selns.push_back(seln("nlna-ht0-met200-mt0-nj0-nb0", /*nleps_*/ -1,  intmax, /*sig_pt*/ 20.,/*veto_pt*/ 10., /*ht*/ 0., dblmax, /*mj*/ 0.,  dblmax, /*met*/ 200.,  dblmax, /*mt*/ 0,  dblmax, /*njets*/0,  intmax, /*nbl*/ 0,  intmax));
  // selns.push_back(seln("nlna-ht0-met0to100-mt0-nj0-nb0", /*nleps_*/ -1,  intmax, /*sig_pt*/ 20.,/*veto_pt*/ 10., /*ht*/ 0., dblmax, /*mj*/ 0.,  dblmax, /*met*/ 0.,  100., /*mt*/ 0,  dblmax, /*njets*/0,  intmax, /*nbl*/ 0,  intmax));
  // selns.push_back(seln("nlna-ht500-met0-mt0-nj0-nb0", /*nleps_*/ -1,  intmax, /*sig_pt*/ 20.,/*veto_pt*/ 10., /*ht*/ 500., dblmax, /*mj*/ 0.,  dblmax, /*met*/ 0., dblmax, /*mt*/ 0,  dblmax, /*njets*/0,  intmax, /*nbl*/ 0,  intmax));
  
  //selns.push_back(seln("nlna-ht500-met200-mt0-nj6-nb0-mj300to450", /*nleps_*/ -1,  intmax, /*sig_pt*/ 20.,/*veto_pt*/ 10., /*ht*/ 500., dblmax, /*mj*/ 300.,  450., /*met*/ 200.,  dblmax, /*mt*/ 0,  dblmax, /*njets*/6,  intmax, /*nbl*/ 0,  intmax));
  //selns.push_back(seln("nlna-ht500-met200-mt0-nj6-nb0-mj450to600", /*nleps_*/ -1,  intmax, /*sig_pt*/ 20.,/*veto_pt*/ 10., /*ht*/ 500., dblmax, /*mj*/ 450.,  600., /*met*/ 200.,  dblmax, /*mt*/ 0,  dblmax, /*njets*/6,  intmax, /*nbl*/ 0,  intmax));
  //selns.push_back(seln("nlna-ht500-met200-mt0-nj6-nb0-mj600toInf", /*nleps_*/ -1,  intmax, /*sig_pt*/ 20.,/*veto_pt*/ 10., /*ht*/ 500., dblmax, /*mj*/ 600., dblmax, /*met*/ 200.,  dblmax, /*mt*/ 0,  dblmax, /*njets*/6,  intmax, /*nbl*/ 0,  intmax));
 
  // selns.push_back(seln("nl1-ht500-met200-mt150-nj6-nb0-mj300to400", /*nleps_*/ 1,  intmax, /*sig_pt*/ 20.,/*veto_pt*/ 10., /*ht*/ 500., dblmax, /*mj*/ 300.,  400., /*met*/ 200.,  dblmax, /*mt*/ 150,  dblmax, /*njets*/6,  intmax, /*nbl*/ 0,  intmax));
  // selns.push_back(seln("nl1-ht500-met200-mt150-nj6-nb0-mj400to500", /*nleps_*/ 1,  intmax, /*sig_pt*/ 20.,/*veto_pt*/ 10., /*ht*/ 500., dblmax, /*mj*/ 400.,  500., /*met*/ 200.,  dblmax, /*mt*/ 150,  dblmax, /*njets*/6,  intmax, /*nbl*/ 0,  intmax));
  // selns.push_back(seln("nl1-ht500-met200-mt150-nj6-nb0-mj500to600", /*nleps_*/ 1,  intmax, /*sig_pt*/ 20.,/*veto_pt*/ 10., /*ht*/ 500., dblmax, /*mj*/ 500.,  600., /*met*/ 200.,  dblmax, /*mt*/ 150,  dblmax, /*njets*/6,  intmax, /*nbl*/ 0,  intmax));
  // selns.push_back(seln("nl1-ht500-met200-mt150-nj6-nb0-mj600toInf", /*nleps_*/ 1,  intmax, /*sig_pt*/ 20.,/*veto_pt*/ 10., /*ht*/ 500., dblmax, /*mj*/ 600.,  dblmax, /*met*/ 200.,  dblmax, /*mt*/ 150,  dblmax, /*njets*/6,  intmax, /*nbl*/ 0,  intmax));
  msgsvc(msg::dbg,"Defined selections...");

  //---------------- DEFINE VARIABLES FOR ALL PLOTS ----------------
  map<TString, var> vars_map;
  vars_map["mj"] = var("mj", "M_{J} [GeV]", 30, 0., 1500.);
  vars_map["gen_mj"] = var("gen_mj", "gen-level M_{J} [GeV]", 30, 0., 1500.);
  vars_map["ht"] = var("ht", "H_{T} [GeV]", 35, 500., 4000.);
  vars_map["ht_isr_me"] = var("ht_isr_me", "H_{T}(ISR ME) [GeV]", 15, 0., 1500.);
  vars_map["ht_fsr"] = var("ht_fsr", "H_{T}(FSR PS) [GeV]", 20, 0., 500.);
  // vars_map["ht_fsr_topdiff"] = var("ht_fsr", "H_{T}(FSR PS) [GeV]", 20, 0., 500.);
  vars_map["ht_top_daughters"] = var("ht_top_daughters", "H_{T}(tops daughters) [GeV]", 20, 0., 2000.);
  vars_map["ht_part"] = var("ht_part", "H_{T}(partons) [GeV]", 35, 500., 4000.);
  vars_map["njets"] = var("njets", "jet multiplicity", 20, 0.5, 20.5);
  vars_map["nbl"] = var("nbl", "b-jet multiplicity", 7, -0.5, 6.5);
  vars_map["fjm1"] = var("fjm1", "m(J_{1}) [GeV]", 10, 0., 500.);
  vars_map["fjm1_nconst"] = var("fjm1_nconst", "m_{j1} const. multiplicity", 5, 0., 5.);
  vars_map["fjm2"] = var("fjm2", "m_{j2} [GeV]", 10, 0., 500.);
  vars_map["fjm2_nconst"] = var("fjm2_nconst", "m_{j2} const. multiplicity", 5, 0., 5.);
  vars_map["fjpt1"] = var("fjpt1", "lead-m large-R jet p_{T} [GeV]", 30, 0., 1500.);
  vars_map["nfjets"] = var("nfjets", "large-R jet multiplicity", 6, 0.5, 6.5);
  vars_map["met"] = var("met", "MET [GeV]", 16, 0., 800.);
  vars_map["met_me"] = var("met_me", "MET (ME) [GeV]", 16, 0., 800.);
  vars_map["mt"] = var("mt", "m_{T} [GeV]", 20, 0., 1000.);
  vars_map["isr"] = var("isr", "ISR [GeV]", 10, 0., 1000.);
  // vars_map["pt_isr"] = var("pt_isr", "ISR p_{T}[GeV]", 100, 0., 1000.);
  // vars_map["pt_gg"] = var("pt_gg", "gluino system p_{T}[GeV]", 100, 0., 1000.);
  // vars_map["lead_isr_pt"] = var("lead_isr_pt", "leading ISR p_{T} [GeV]", 20, 0., 1000.);
  // vars_map["sublead_isr_pt"] = var("sublead_isr_pt", "subleading ISR p_{T} [GeV]", 20, 0., 1000.);
  vars_map["nisr_me"] = var("nisr_me", "ISR ME parton multiplicity", 5, -0.5, 4.5);
  // vars_map["npart"] = var("npart", "parton multiplicity", 20, 0.5, 20.5);
  // vars_map["ntop_daughters"] = var("ntop_daughters", "top daughter multiplicity", 7, 0.5, 7.5);
  vars_map["ptt"] = var("ptt", "p_{T}(t#bar{t}), ME [GeV]", 100, 0., 1000.);
  vars_map["dphi_tt"] = var("dphi_tt", "#Delta#phi^{gen}(t,#bar{t})", 6, 0., 3.14);
  vars_map["dphi_fjm1_fjm2"] = var("dphi_fjm1_fjm2", "#Delta#phi^{gen}(m_{j1},m_{j2})", 10, 0., 3.14);
  // vars_map["lead_pt_top"] = var("lead_pt_top", "leading top p^{gen}_{T} [GeV]", 30, 0., 1500.);
  // vars_map["sublead_pt_top"] = var("sublead_pt_top", "sub-leading top p^{gen}_{T} [GeV]", 30, 0., 1500.);
  msgsvc(msg::dbg,"Defined variables...");

  //---------------- DEFINE VARIABLE PAIRS 2D HISTOGRAMS ----------------
  vector<pair<TString,TString> > var_pairs; //use only variables defined in vars_map
  // var_pairs.push_back(make_pair("ht","mj"));
  // var_pairs.push_back(make_pair("njets","mj"));
  // var_pairs.push_back(make_pair("njets","ht"));
  // var_pairs.push_back(make_pair("nbl","mj"));
  // var_pairs.push_back(make_pair("met","mj"));
  // var_pairs.push_back(make_pair("met","ht"));
  // var_pairs.push_back(make_pair("met","njets"));
  // var_pairs.push_back(make_pair("lead_pt_top", "ptt"));
  // var_pairs.push_back(make_pair("lead_pt_top", "fjpt1"));
  // var_pairs.push_back(make_pair("ptt","mj"));
  // var_pairs.push_back(make_pair("lead_pt_top","mj"));
  // var_pairs.push_back(make_pair("dphi_tt","ptt"));
  // var_pairs.push_back(make_pair("dphi_tt","mj"));

  // var_pairs.push_back(make_pair("pt_gg","pt_isr"));
  // var_pairs.push_back(make_pair("ptt","pt_isr"));
  // var_pairs.push_back(make_pair("fjm1","ht_isr_me"));
  // var_pairs.push_back(make_pair("fjm1","nisr_me"));

  // var_pairs.push_back(make_pair("nisr_me","mj"));
  // var_pairs.push_back(make_pair("nisr_me","dphi_tt"));
  // var_pairs.push_back(make_pair("nisr_me","ht_isr_me"));

  // var_pairs.push_back(make_pair("dphi_tt","ht_isr_me"));
  // var_pairs.push_back(make_pair("dphi_tt","nisr_me"));
  // var_pairs.push_back(make_pair("ht_isr_me","mj"));

  var_pairs.push_back(make_pair("met","met_me"));
  var_pairs.push_back(make_pair("mj","gen_mj"));

  msgsvc(msg::dbg,"Defined correlations...");

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

  msgsvc(msg::dbg,"Defined 3D plots...");

  //---------------- LOOP OVER SAMPLES ----------------
  for (vector<sample>::iterator isamp = samples.begin(); isamp != samples.end(); isamp++) {

    if (isamp->name!="ttbar") split_mctype = false;//assuming ttbar is always first in the vector...
    msgsvc(msg::dbg, "Running over sample: " + isamp->name);

    map<TString, h1d> h1d_map;
    map<TString, h2d> h2d_map;
    map<TString, h3d> h3d_map;

    TString outdir = "./ntup/mj_plots_"+isamp->name+".root";
    TFile fout(outdir,"RECREATE");
    small_tree_quick tree((basedir+isamp->filestr).Data());

    //---------------- CREATE 1D HISTOGRAMS FOR ALL VARIABLES ----------------
    for (map<TString, var>::iterator ivar = vars_map.begin(); ivar!=vars_map.end(); ivar++){
      for (vector<seln>::iterator iseln = selns.begin(); iseln!=selns.end(); iseln++){
        if (split_mctype){
          for (vector<TString>::iterator itype = mctypes.begin(); itype!=mctypes.end(); itype++){
            TString hname = iseln->name+"_"+ivar->second.name+"_"+isamp->name + "_" + *itype;
            msgsvc(msg::dbg, "Creating histogram: " + hname);
            h1d_map[hname] = h1d(ivar->second,hname);
          }
        } else {
          TString hname = iseln->name+"_"+ivar->second.name+"_"+isamp->name;
          msgsvc(msg::dbg, "Creating histogram: " + hname);
          h1d_map[hname] = h1d(ivar->second,hname);
        }
      }
    }

    //---------------- CREATE 2D HISTOGRAMS FOR ALL VARIABLE PAIRS ----------------
    for (vector<pair<TString,TString> >::iterator ipair=var_pairs.begin(); ipair!=var_pairs.end(); ipair++){
      for (vector<seln>::iterator iseln = selns.begin(); iseln!=selns.end(); iseln++){
        if (split_mctype){
          for (vector<TString>::iterator itype = mctypes.begin(); itype!=mctypes.end(); itype++){
            TString hname = iseln->name+"_"+ipair->first+"_"+ipair->second+"_"+isamp->name + "_" + *itype;
            msgsvc(msg::dbg, "Creating histogram: " + hname);
            h2d_map[hname] = h2d(vars_map[ipair->first], vars_map[ipair->second], hname);
          }
        } else {
          TString hname = iseln->name+"_"+ipair->first+"_"+ipair->second+"_"+isamp->name;
          msgsvc(msg::dbg, "Creating histogram: " + hname);
          h2d_map[hname] = h2d(vars_map[ipair->first], vars_map[ipair->second], hname);
        }
      }
    }

    //---------------- CREATE 3D HISTOGRAMS FOR ALL VARIABLE SETS ----------------
    for (vector<vector<TString> >::iterator iset=var_sets.begin(); iset!=var_sets.end(); iset++){
      for (vector<seln>::iterator iseln = selns.begin(); iseln!=selns.end(); iseln++){
        if (split_mctype){
          for (vector<TString>::iterator itype = mctypes.begin(); itype!=mctypes.end(); itype++){
            TString hname = iseln->name+"_"+iset->at(0)+"_"+iset->at(1)+"_"+iset->at(2)+"_"+isamp->name + "_" + *itype;
            msgsvc(msg::dbg, "Creating histogram: " + hname);
            h3d_map[hname] = h3d(vars_map[iset->at(0)], vars_map[iset->at(1)], vars_map[iset->at(2)], hname);
          }
        } else {
          TString hname = iseln->name+"_"+iset->at(0)+"_"+iset->at(1)+"_"+iset->at(2)+"_"+isamp->name;
          msgsvc(msg::dbg, "Creating histogram: " + hname);
          h3d_map[hname] = h3d(vars_map[iset->at(0)], vars_map[iset->at(1)], vars_map[iset->at(2)], hname);
        }
      }
    }

    //---------------- LOOP OVER EVENTS ----------------
    const size_t nent = tree.GetEntries();
    msgsvc(msg::info, TString::Format("Number of events to run over: %u", unsigned(nent)));
    for (size_t ientry(0); ientry<nent; ientry++){
      if (ientry%100000==0) cout<<"INFO:: Processed events: "<<ientry<<endl;
      tree.GetEntry(ientry);
      msgsvc(msg::dbg, TString::Format("Event number: %u", unsigned(tree.event())));
      unsigned clean_mctype = tree.mc_type()&0x0FFF; 
      TString mctype = getMCType(clean_mctype);

      // if ((mctype&&0x0FFF)!=0x101 && (mctype&&0x0FFF)!=0x112) continue;

      double weight = 1.;//tree.weight()*5;
      // events that have the corresponding event display
      // if (tree.event()!= 106827864 && 
      //     tree.event()!=105178777 && 
      //     tree.event()!=73140398 && 
      //     tree.event()!=67007130 && 
      //     tree.event()!=59229595 && 
      //     tree.event()!=45145334 && 
      //     tree.event()!=23202067) continue;

      //mismeasured met
      // if (tree.event()!=85869132 &&
      //     tree.event()!=93026692 &&
      //     tree.event()!=15423267 &&
      //     tree.event()!=58977266 &&
      //     tree.event()!=9011482 &&
      //     tree.event()!=78418880) continue;

      // msglvl = msg::truth;
      // if (ientry!=90000) continue;

      //------------- PRESELECTION -------------------
      if (tree.ht()<500. || tree.met()<200.) continue;

      map<TString, double> flt_val;
      // for the sake of filling histos all are floats...
      flt_val["ht"] = tree.ht();
      flt_val["njets"] = double(tree.njets()); 
      flt_val["nbl"] = double(tree.nbl()); 
      flt_val["mj"] = tree.mj();
      flt_val["met"] = tree.met(); 
      flt_val["mt"] = tree.mt(); 
      flt_val["nfjets"] = tree.nfjets(); 

      // ------------ CALCULATE CUSTOM VARIABLES --------------
      msgsvc(msg::dbg, "Ordering fat jets ");
      // ------ MJ with massive fat jets and leading mass jet vars ----------
      flt_val["fjm1"] = 0.; 
      flt_val["fjpt1"] = 0.; 
      if (tree.fjets_m().size()>0){
        vector<TLorentzVector> fjs;
        for (size_t i(0); i<tree.fjets_m().size(); i++){
          TLorentzVector ifj;
          ifj.SetPtEtaPhiM(tree.fjets_pt()[i],tree.fjets_eta()[i],tree.fjets_phi()[i],tree.fjets_m()[i]);; 
          fjs.push_back(ifj);        
        }
        flt_val["fjpt1"] = fjs[0].Pt();
        sort(fjs.begin(), fjs.end(), compMass);
        flt_val["fjm1"] = fjs[0].M();
      }
      
      msgsvc(msg::dbg, "Looking at the generator level ");
      // ------ Generator level ----------
      vector<TVector3> tops_me; 
      vector<unsigned> tops_me_ind; 
      vector<TVector3> gluinos; 
      vector<TVector3> tops_ps; 
      vector<unsigned> tops_ps_ind; 
      vector<TVector3> top_daughters; 
      set<unsigned> top_daughters_ind; 
      vector<TVector3> top_bs; 
      vector<TVector3> isr_me; 
      vector<TVector3> fsr; 
      TVector3 met_me; 
      if (msglvl == msg::truth) cout<<"==================== Event: "<<tree.event()<<" ======================="<<endl;
      // int ancestor_id = intmax;
      // for (size_t i(0); i<tree.mc_id().size(); i++){
      //   if (isamp->name=="ttbar"){
      //     if (abs(tree.mc_id()[i])==6 && tree.mc_status()[i]==22) ancestor_id = tree.mc_mom()[i];
      //   } 
      // }

      unsigned nleps = 0;
      unsigned ntaus = 0;
      unsigned nvs = 0; // # neutrinos
      for (size_t i(0); i<tree.mc_id().size(); i++){
        //don't care about the neutralino
        if (abs(tree.mc_id()[i])==1000022) continue;

        //clean out some garbage ...
        if (tree.mc_status()[i]==52 || tree.mc_status()[i]==0 || 
            (tree.mc_status()[i]>1 && tree.mc_status()[i]<10) || 
            (tree.mc_status()[i]==1 && abs(tree.mc_id()[i])>22) ||
            (tree.mc_status()[i]>70 && tree.mc_status()[i]<74)) {
          if (abs(tree.mc_id()[i])!=15) continue; 
        }

        if (msglvl == msg::truth) cout<<setw(2)<<i<<": Id  "<<setw(10)<<tree.mc_id()[i]<<
              "  Sta "<<setw(10)<<tree.mc_status()[i]<<
              "  Mom "<<setw(10)<<tree.mc_mom()[i]<<
              "  Eta "<<setw(10)<<tree.mc_eta()[i]<<
              "  Phi "<<setw(10)<<tree.mc_phi()[i]<<
              "  Pt "<<setw(10)<<tree.mc_pt()[i];
        // if (tree.mc_pt()[i]<30. || fabs(tree.mc_eta()[i])>5.) {
              // if (abs(tree.mc_id()[i])!=12 && abs(tree.mc_id()[i])!=14 && abs(tree.mc_id()[i])!=16) {
        //   if (msglvl == msg::truth) cout<<" Outside acceptance"<<endl;
        //   continue;
              // }
        // }

        // tops from the matrix element
        if (abs(tree.mc_id()[i])==1000021) {
          TVector3 iglu; iglu.SetPtEtaPhi(tree.mc_pt()[i], tree.mc_eta()[i], tree.mc_phi()[i]);
          gluinos.push_back(iglu);
          if (msglvl == msg::truth) cout<<" Gluino"<<endl;

        } else if (abs(tree.mc_id()[i])==6 && tree.mc_status()[i]==22) {
          TVector3 itop_me; itop_me.SetPtEtaPhi(tree.mc_pt()[i], tree.mc_eta()[i], tree.mc_phi()[i]);
          tops_me.push_back(itop_me);
          tops_me_ind.push_back(i);
          if (msglvl == msg::truth) cout<<" ME top"<<endl;

        // tops that decay to W and b
        } else if (abs(tree.mc_id()[i])==6 && tree.mc_status()[i]!=22) {
          TVector3 itop_ps; itop_ps.SetPtEtaPhi(tree.mc_pt()[i], tree.mc_eta()[i], tree.mc_phi()[i]);
          tops_ps.push_back(itop_ps);
          tops_ps_ind.push_back(i);
          if (msglvl == msg::truth) cout<<" PS top"<<endl;

        // b's from top
        } else if (abs(tree.mc_id()[i])==5 && abs(tree.mc_id()[tree.mc_mom()[i]])==6 && tree.mc_status()[i]==23){
          TVector3 itop_child; itop_child.SetPtEtaPhi(tree.mc_pt()[i], tree.mc_eta()[i], tree.mc_phi()[i]);
          top_daughters.push_back(itop_child);
          top_daughters_ind.insert(i);
          top_bs.push_back(itop_child);
          if (msglvl == msg::truth) cout<<" Top child"<<endl;
        
        // W from top
        } else if (abs(tree.mc_id()[i])==24 && abs(tree.mc_id()[tree.mc_mom()[i]])==6  && tree.mc_status()[i]==22){
          if (msglvl == msg::truth) cout<<" W boson"<<endl;
          continue;
          
        // categorize non-tau W children
        } else if (abs(tree.mc_id()[tree.mc_mom()[i]])==24 && (tree.mc_status()[i]==23 || tree.mc_status()[i]==1)){
          //remove duplicates caused by taking particles from both mc_doc and mc_final
          bool duplicate = false;
          for (size_t j(0); j<i; j++){
            if (tree.mc_id()[i] == tree.mc_id()[j] && 
                fabs(tree.mc_pt()[i]-tree.mc_pt()[j])<0.0001  && 
                fabs(tree.mc_eta()[i]-tree.mc_eta()[j])<0.0001  && 
                fabs(tree.mc_phi()[i]-tree.mc_phi()[j])<0.0001) {
              duplicate = true;
              break;
            }
          }
          if (duplicate) {
            if (msglvl == msg::truth) cout<<" Duplicate"<<endl;
            continue;
          }

          if (abs(tree.mc_id()[i])==12 || abs(tree.mc_id()[i])==14 || abs(tree.mc_id()[i])==16) {
            if (msglvl == msg::truth) cout<<" Neutrino"<<endl;
            if (nvs==0) {
              met_me.SetPtEtaPhi(tree.mc_pt()[i], tree.mc_eta()[i], tree.mc_phi()[i]);
            } else {
              TVector3 inu; inu.SetPtEtaPhi(tree.mc_pt()[i], tree.mc_eta()[i], tree.mc_phi()[i]);
              met_me = met_me + inu;
            }
            nvs++;
            continue;
          } else if (abs(tree.mc_id()[i])==11 || abs(tree.mc_id()[i])==13){
            nleps++;
          } else if (abs(tree.mc_id()[i])==15){
            ntaus++;
          }
          TVector3 itop_child; itop_child.SetPtEtaPhi(tree.mc_pt()[i], tree.mc_eta()[i], tree.mc_phi()[i]);
          top_daughters.push_back(itop_child);
          top_daughters_ind.insert(i);
          if (msglvl == msg::truth) cout<<" Top child"<<endl;

        //classify taus
        } else if (abs(tree.mc_id()[i])==15 && abs(tree.mc_id()[tree.mc_mom()[i]])==24){
          ntaus++;
          TVector3 itop_child; itop_child.SetPtEtaPhi(tree.mc_pt()[i], tree.mc_eta()[i], tree.mc_phi()[i]);
          top_daughters.push_back(itop_child);
          top_daughters_ind.insert(i);
          if (msglvl == msg::truth) cout<<" Top child"<<endl;

        //classify tau children
        } else if (abs(tree.mc_id()[tree.mc_mom()[i]])==15){
          bool duplicate = false;
          for (size_t j(0); j<i; j++){
            if (tree.mc_id()[i] == tree.mc_id()[j] && 
                fabs(tree.mc_pt()[i]-tree.mc_pt()[j])<0.0001  && 
                fabs(tree.mc_eta()[i]-tree.mc_eta()[j])<0.0001  && 
                fabs(tree.mc_phi()[i]-tree.mc_phi()[j])<0.0001) {
              duplicate = true;
              break;
            }
          }
          if (duplicate) {
            if (msglvl == msg::truth) cout<<" Duplicate"<<endl;
            continue;
          }

          if (abs(tree.mc_id()[i])==12 || abs(tree.mc_id()[i])==14 || abs(tree.mc_id()[i])==16) {
            if (msglvl == msg::truth) cout<<" Neutrino from tau"<<endl;
            if (nvs==0) {
              met_me.SetPtEtaPhi(tree.mc_pt()[i], tree.mc_eta()[i], tree.mc_phi()[i]);
            } else {
              TVector3 inu; inu.SetPtEtaPhi(tree.mc_pt()[i], tree.mc_eta()[i], tree.mc_phi()[i]);
              met_me = met_me + inu;
            }
            nvs++;
            continue;
          }
          if (msglvl == msg::truth) cout<<" Tau child - don't count"<<endl;

        // fsr? 
        } else if(abs(tree.mc_id()[tree.mc_mom()[i]])==6 && abs(tree.mc_id()[i])!=24 && abs(tree.mc_id()[i])!=5) {
          TVector3 ifsr; ifsr.SetPtEtaPhi(tree.mc_pt()[i], tree.mc_eta()[i], tree.mc_phi()[i]);
          fsr.push_back(ifsr);
          if (msglvl == msg::truth) cout<<" FSR"<<endl;

        // ISR from ME
        } else if (tree.mc_status()[i]==23) {
          TVector3 iisr_me; iisr_me.SetPtEtaPhi(tree.mc_pt()[i], tree.mc_eta()[i], tree.mc_phi()[i]);
          isr_me.push_back(iisr_me);
          if (msglvl == msg::truth) cout<<" ME ISR"<<endl;

        // odd left-overs...? 
        } else {
          //ignore additional photon FSR for now...
          if (abs(tree.mc_id()[i])==22) {
            if (msglvl == msg::truth) cout<<" Duplicate ... "<<endl;
            continue; 
          //more duplicates  
          } else if (tree.mc_status()[i]==1 && tree.mc_id()[tree.mc_mom()[i]]==tree.mc_id()[i]) {
            if (msglvl == msg::truth) cout<<" Duplicate ... "<<endl;
            continue; 
          //mostly garbage attached to the proton, this is not ME ISR or top so it is from the parton shower
          } else if (tree.mc_mom()[i]>5000) {
            if (msglvl == msg::truth) cout<<" Garbage ... "<<endl;
            continue; 
          }
          cout<<" **** ??? **** "<<tree.event()<<endl;
        }
      }  

      bool found_dupl = false;
      for (unsigned itop(0);itop<tops_ps.size();itop++){
        for (unsigned jtop(0);jtop<tops_ps.size();jtop++){
          if (tree.mc_mom()[tops_ps_ind[itop]]==tree.mc_mom()[tops_ps_ind[jtop]]){
            if (tree.mc_pt()[tops_ps_ind[itop]]<tree.mc_pt()[tops_ps_ind[jtop]]) {
              tops_ps.erase(tops_ps.begin()+itop);
              tops_ps_ind.erase(tops_ps_ind.begin()+itop);
              found_dupl = true;
              break;
            } 
          }
          if (found_dupl) break;
        }
      }

      // TVector3 diff1, diff2;
      // if (tops_me.size()==2 && tops_ps.size()==2){
      //   if (tree.mc_mom()[tops_ps_ind[0]]==tops_me_ind[0]){
      //     if (tree.mc_mom()[tops_ps_ind[1]]==tops_me_ind[1]){
      //       diff1 = tops_me[0] - tops_ps[0];
      //       diff2 = tops_me[1] - tops_ps[1];
      //     } else {
      //       cout<<"ERROR 1"<<endl;
      //     }
      //   } else if (tree.mc_mom()[tops_ps_ind[0]]==tops_me_ind[1]){
      //     if (tree.mc_mom()[tops_ps_ind[1]]==tops_me_ind[0]){
      //       diff1 = tops_me[0] - tops_ps[1];
      //       diff2 = tops_me[1] - tops_ps[0];
      //     } else {
      //       cout<<"ERROR 2"<<endl;
      //     }
      //   } else {
      //     cout<<"ERROR"<<endl;
      //   }
      // } else {
      //   if (msglvl == msg::info) cout<<"Found "<<tops_me.size()<<" ME tops and "<<tops_ps.size()<<" PS tops"<<endl;
      //   for (unsigned i(0); i< tops_ps.size(); i++){
      //     cout<<"Top with mom "<<tree.mc_mom()[tops_ps_ind[i]]<< " has pt "<<tree.mc_pt()[tops_ps_ind[i]]<<endl;
      //   }
      // }
      // flt_val["ht_fsr_topdiff"] = diff1.Pt() + diff2.Pt();

      // flt_val["ntop_daughters"] = top_daughters.size();
      flt_val["ntops_me"] = tops_me.size();
      flt_val["ntops_ps"] = tops_ps.size();
      flt_val["nisr_me"] = isr_me.size();
      flt_val["npart"] = isr_me.size() + top_daughters.size();

      flt_val["lead_isr_pt"] = -1.;
      flt_val["sublead_isr_pt"] = -1.;
      if (isr_me.size()>0) flt_val["lead_isr_pt"] = isr_me[0].Pt();
      if (isr_me.size()>1) flt_val["sublead_isr_pt"] = isr_me[1].Pt();

      flt_val["met_me"] = met_me.Pt();

      flt_val["ht_isr_me"] = getHT(isr_me);
      flt_val["ht_fsr"] = getHT(fsr);
      flt_val["ht_top_daughters"] = getHT(top_daughters);

      size_t ntops_me = tops_me.size();
      // sort(tops_me.begin(), tops_me.end(), compPT);
      // msgsvc(msg::dbg, "Record top pTs ");
      // flt_val["lead_pt_top"] = -1.*dblmax;
      // flt_val["sublead_pt_top"] = -1.*dblmax;
      // if (ntops_me>0) flt_val["lead_pt_top"] = tops_me[0].Pt();
      // if (ntops_me>1) flt_val["sublead_pt_top"] = tops_me[1].Pt();

      // ------ Generator level pt of the ttbar/gg system and the isr - these should balance ----------
      if (ntops_me!=2){
        if (ntops_me!=2 && isamp->name=="ttbar") msgsvc(msg::err, "Bad truth record, number of tops found = " + TString::Format("%i",int(ntops_me)));
        flt_val["ptt"] = dblmax; 
        flt_val["dphi_tt"] = dblmax;
      } else {
        TVector3 tt_me = tops_me[0]+tops_me[1];
        flt_val["ptt"] = tt_me.Pt();
        flt_val["dphi_tt"] = tops_me[0].DeltaPhi(tops_me[1]);
      }

      // flt_val["pt_gg"] = dblmax;
      // if (isamp->name=="T1tttt1200" || isamp->name=="T1tttt1500"){
      //   if (gluinos.size()==2){
      //     TVector3 glu_pair = gluinos[0]+gluinos[1];
      //     flt_val["pt_gg"] = glu_pair.Pt();
      //   } else {
      //     msgsvc(msg::err, "Bad truth record, number of gluinos found = " + TString::Format("%i",int(gluinos.size())));
      //   }
      // } 

      // flt_val["pt_isr"] = dblmax;
      // TVector3 isr_system;
      // for (size_t iisr(0); iisr<isr_me.size(); iisr++){
      //   if (iisr==0) {
      //     isr_system.SetPtEtaPhi(isr_me[iisr].Pt(),isr_me[iisr].Eta(), isr_me[iisr].Phi());
      //   } else {
      //     TVector3 v_iisr; v_iisr.SetPtEtaPhi(isr_me[iisr].Pt(),isr_me[iisr].Eta(), isr_me[iisr].Phi());
      //     isr_system += v_iisr;
      //   }
      // }
      // flt_val["pt_isr"] = isr_system.Pt();

      //------------------------- Try to deduce the FSR from looking at the diff of the top ME and top PS

      msgsvc(msg::dbg, "Beginning to fill histograms");

      for (vector<seln>::iterator iseln = selns.begin(); iseln!=selns.end(); iseln++){
        msgsvc(msg::dbg, "Checking selection "+ iseln->name);      
        if (!passSelection(tree, *iseln)) continue;
        msgsvc(msg::dbg, "Selection "+iseln->name + " passed");

        for (map<TString, var>::iterator ivar = vars_map.begin(); ivar!=vars_map.end(); ivar++){
          msgsvc(msg::dbg, "Filling "+ivar->first);
          TString hname = iseln->name+"_"+ivar->second.name+"_"+isamp->name;
          if (split_mctype) hname = iseln->name+"_"+ivar->second.name+"_"+isamp->name + "_" + mctype;
          h1d_map[hname].fill(flt_val[ivar->second.name], weight);
        }

        for (vector<pair<TString,TString> >::iterator ipair=var_pairs.begin(); ipair!=var_pairs.end(); ipair++){
          TString hname = iseln->name+"_"+ipair->first+"_"+ipair->second+"_"+isamp->name;
          if (split_mctype) hname = iseln->name+"_"+ipair->first+"_"+ipair->second+"_"+isamp->name + "_" + mctype;
          h2d_map[hname].fill(flt_val[ipair->first], flt_val[ipair->second], weight);
        }

        // for (vector<vector<TString> >::iterator iset=var_sets.begin(); iset!=var_sets.end(); iset++){
        //   TString hname = iseln->name+"_"+iset->at(0)+"_"+iset->at(1)+"_"+iset->at(2)+"_"+isamp->name;
        //   if (split_mctype) hname = iseln->name+"_"+iset->at(0)+"_"+iset->at(1)+"_"+iset->at(2)+"_"+isamp->name + "_" + mctype;
        //   h3d_map[hname].fill(flt_val[iset->at(0)], flt_val[iset->at(1)], flt_val[iset->at(2)], weight);
        // }

      }    
    }

    fout.Write();
    fout.Close();
  }

  return 0;
}

bool compMass(const TLorentzVector &tlv1, const TLorentzVector &tlv2){
  return (tlv1.M()>tlv2.M());
}

bool compPT(const TVector3 &v1, const TVector3 &v2){
  return (v1.Pt()>v2.Pt());
}

double getHT(vector<TVector3> &v){
  double ht = 0.;
  for (size_t i(0); i<v.size(); i++) ht +=v[i].Pt();
  return ht;
}

bool passSelection(const small_tree &tree, const seln &iseln){

  size_t nveto_mus(0), nveto_els(0);
  vector<int> sigels_index = vector<int>(0);
  vector<int> sigmus_index = vector<int>(0);
  msgsvc(msg::dbg, "Check selection");

  if (iseln.nleps>=0){ //if we care about reco leptons
    //---------- ELECTRONS ----------------      
    for (size_t iel=0; iel<tree.els_pt().size(); iel++){
      if (!tree.els_ispf()[iel]) continue;
      if (tree.els_pt()[iel] < iseln.veto_lep_pt_min) continue;
      if (tree.els_miniso_tr10()[iel] < 0.1) nveto_els++;


      if (!(tree.els_sigid()[iel])) continue;
      if (tree.els_pt()[iel] < iseln.sig_lep_pt_min) continue;
      if (tree.els_miniso_tr10()[iel] < 0.1) sigels_index.push_back(iel);
    }
    //---------- MUONS ----------------
    for (size_t imu=0; imu<tree.mus_pt().size(); imu++){
      if (tree.mus_pt()[imu] < iseln.veto_lep_pt_min) continue;
      if (tree.mus_miniso_tr10()[imu] < 0.2) nveto_mus++;

      if (!(tree.mus_sigid()[imu])) continue;
      if (tree.mus_pt()[imu] < iseln.sig_lep_pt_min) continue;
      if (tree.mus_miniso_tr10()[imu] < 0.2) sigmus_index.push_back(imu);
    }

    if ((sigels_index.size()+sigmus_index.size())!=unsigned(iseln.nleps) || (nveto_els+nveto_mus)!=unsigned(iseln.nleps)) return false;
  }

  if (iseln.nleps==1) {
    msgsvc(msg::dbg, "Calculating mT");
    // double mt = 0;
    // if (sigels_index.size()>0) mt = sqrt(2*tree.els_pt()[sigels_index[0]]*tree.met()*(1-cos(tree.met_phi()-tree.els_phi()[sigels_index[0]])));
    // else if (sigmus_index.size()>0) mt = sqrt(2*tree.mus_pt()[sigmus_index[0]]*tree.met()*(1-cos(tree.met_phi()-tree.mus_phi()[sigmus_index[0]])));
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

TString getMCType(unsigned mc_type){
  // if (mc_type==0x200 || mc_type==0x211 || mc_type==0x222) return "ll";
  // else if (mc_type==0x101 || mc_type==0x112) return "lt";
  // else if (mc_type==0x100 || mc_type==0x111) return "lh";
  // else if (mc_type==0x000 || mc_type==0x001 || mc_type==0x002) return "hh";
  // else if (mc_type==0x000) return "hh";
  unsigned nleps = (mc_type&0xF00)>>8;
  nleps += (mc_type&0x00F);
  nleps -= (mc_type&0x0F0)>>4;

  if (nleps==2) return "nl2";
  else if (nleps==1) return "nl1";
  else if (nleps==0) return "nl0";

  return "";
}

void msgsvc(const unsigned &lvl, const TString &mymsg) {
  if (lvl>=msglvl) {
    if(lvl == msg::dbg) cout<<"DEBUG:: "<<mymsg<<endl;
    else if(lvl == msg::info) cout<<"INFO:: "<<mymsg<<endl;
    else if(lvl == msg::warn) cout<<"WARNING:: "<<mymsg<<endl;
    else if(lvl == msg::err) cout<<"ERROR:: "<<mymsg<<endl;

  }
}
