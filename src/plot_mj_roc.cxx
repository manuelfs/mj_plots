// plot_roc: Macro that plots ROC curves

#include <stdexcept>
#include <iostream>

#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TDirectory.h"
#include "TMarker.h"

#include "styles.hpp"
#include "utilities.hpp"
#include "plot_mj_roc.hpp"

using namespace std;

int main(){
  styles style("1Dtitle"); style.setDefaultStyle();

  TString folder="/afs/cern.ch/user/m/manuelf/work/ucsb/15-01-30_skim/";

  // NTUPLES
  vector<TString> v_t1t;
  v_t1t.push_back(folder+"*T1tttt*1500_*PU20*");

  vector<TString> v_t1tc;
  v_t1tc.push_back(folder+"*T1tttt*1200_*PU20*");

  vector<TString> v_t1q;
  v_t1q.push_back(folder+"*T1qqqq*1400_*PU20*");

  vector<TString> v_t1qc;
  v_t1qc.push_back(folder+"*T1qqqq*1000_*PU20*");

  vector<TString> v_tt;
  v_tt.push_back(folder+"*TTJet*");

  vector<TString> v_qcd;
  v_qcd.push_back(folder+"*QCD_HT*");
  //v_qcd.push_back(folder+"*QCD_Pt*");

  vector<TString> v_zjets;
  v_zjets.push_back(folder+"*ZJetsToNuNu*");

  ////////////////////////// SAMPLES for the axes /////////////////////////
  vector<sample_class> tt_t1t; 
  tt_t1t.push_back(sample_class("T1tttt(1500,100)", v_t1t));
  tt_t1t.push_back(sample_class("t#bar{t}", v_tt));

  vector<sample_class> qcd_t1t; 
  qcd_t1t.push_back(sample_class("T1tttt(1500,100)", v_t1t));
  qcd_t1t.push_back(sample_class("QCD", v_qcd));

  vector<sample_class> zjets_t1t; 
  zjets_t1t.push_back(sample_class("T1tttt(1500,100)", v_t1t));
  zjets_t1t.push_back(sample_class("Z#rightarrow#nu#nu", v_zjets));

  vector<sample_class> tt_t1tc; 
  tt_t1tc.push_back(sample_class("T1tttt(1200,800)", v_t1tc));
  tt_t1tc.push_back(sample_class("t#bar{t}", v_tt));

  vector<sample_class> qcd_t1tc; 
  qcd_t1tc.push_back(sample_class("T1tttt(1200,800)", v_t1tc));
  qcd_t1tc.push_back(sample_class("QCD", v_qcd));

  vector<sample_class> zjets_t1tc; 
  zjets_t1tc.push_back(sample_class("T1tttt(1200,800)", v_t1tc));
  zjets_t1tc.push_back(sample_class("Z#rightarrow#nu#nu", v_zjets));

  // T1qqqq
  vector<sample_class> tt_t1q; 
  tt_t1q.push_back(sample_class("T1qqqq(1400,100)", v_t1q));
  tt_t1q.push_back(sample_class("t#bar{t}", v_tt));

  vector<sample_class> qcd_t1q; 
  qcd_t1q.push_back(sample_class("T1qqqq(1400,100)", v_t1q));
  qcd_t1q.push_back(sample_class("QCD", v_qcd));

  vector<sample_class> zjets_t1q; 
  zjets_t1q.push_back(sample_class("T1qqqq(1400,100)", v_t1q));
  zjets_t1q.push_back(sample_class("Z#rightarrow#nu#nu", v_zjets));

  vector<sample_class> tt_t1qc; 
  tt_t1qc.push_back(sample_class("T1qqqq(1000,800)", v_t1qc));
  tt_t1qc.push_back(sample_class("t#bar{t}", v_tt));

  vector<sample_class> qcd_t1qc; 
  qcd_t1qc.push_back(sample_class("T1qqqq(1000,800)", v_t1qc));
  qcd_t1qc.push_back(sample_class("QCD", v_qcd));

  vector<sample_class> zjets_t1qc; 
  zjets_t1qc.push_back(sample_class("T1qqqq(1000,800)", v_t1qc));
  zjets_t1qc.push_back(sample_class("Z#rightarrow#nu#nu", v_zjets));



  ///////////////////// VARIABLES for each ROC /////////////////////
  vector<var_class> mj_general;
  mj_general.push_back(var_class("ht",4000,0,"H_{T}^{40}",8,1));
  mj_general.push_back(var_class("met",1500,0,"MET",2,1));
  mj_general.push_back(var_class("mj_30",2000,0,"M_{J}",4,1));

  vector<var_class> mj_cands;
  //mj_cands.push_back(var_class("ht",4000,0,"H_{T}^{40}",2,1));
  mj_cands.push_back(var_class("mj_cands",2200,0,"M_{J} pfcands",2,1));
  mj_cands.push_back(var_class("mj_cands_trim",2200,0,"M_{J} pfcands trimmed",28,1));
  mj_cands.push_back(var_class("mj_30", 2200,0,"M_{J} 30 GeV jets",4,1));

  vector<var_class> mj_sizes;
  //mj_sizes.push_back(var_class("ht",4000,0,"H_{T}^{40}",2,1));
  mj_sizes.push_back(var_class("mj_r08",2200,0,"M_{J} R=0.8",2,1));
  mj_sizes.push_back(var_class("mj_r10",2200,0,"M_{J} R=1.0",28,1));
  mj_sizes.push_back(var_class("mj_30", 2200,0,"M_{J} R=1.2",4,1));
  mj_sizes.push_back(var_class("mj_r14",2200,0,"M_{J} R=1.4",kMagenta+2,1));

  vector<var_class> mj_pt;
  //mj_pt.push_back(var_class("ht",4000,0,"H_{T}^{40}",2,1));
  mj_pt.push_back(var_class("mj_10",2200,0,"M_{J} 10 GeV jets",2,1));
  mj_pt.push_back(var_class("mj_20",2200,0,"M_{J} 20 GeV jets",28,1));
  mj_pt.push_back(var_class("mj_30", 2200,0,"M_{J} 30 GeV jets",4,1));
  mj_pt.push_back(var_class("mj_40",2200,0,"M_{J} 40 GeV jets",kMagenta+2,1));

  vector<var_class> mj_eta;
  //mj_eta.push_back(var_class("ht",4000,0,"H_{T}^{40}",2,1));
  mj_eta.push_back(var_class("mj_eta25",2200,0,"M_{J} |#eta| < 2.5",2,1));
  mj_eta.push_back(var_class("mj_30", 2200,0,"M_{J} |#eta| < 5",4,1));

  vector<var_class> mj_ptfat;
  //mj_ptfat.push_back(var_class("ht",4000,0,"H_{T}^{40}",2,1));
  mj_ptfat.push_back(var_class("Sum$(fjets_30_m*(fjets_30_pt>=30))",2200,0,"M_{J} 30 GeV large-R jets",2,1));
  mj_ptfat.push_back(var_class("mj_30", 2200,0,"M_{J} 50 GeV large-R jets",4,1));
  mj_ptfat.push_back(var_class("Sum$(fjets_30_m*(fjets_30_pt>=70))",2200,0,"M_{J} 70 GeV large-R jets",28,1));
  mj_ptfat.push_back(var_class("Sum$(fjets_30_m*(fjets_30_pt>=100))",2200,0,"M_{J} 100 GeV large-R jets",kMagenta+2,1));

  vector<var_class> mj_lep;
  //mj_lep.push_back(var_class("ht",4000,0,"H_{T}^{40}",2,1));
  mj_lep.push_back(var_class("mj_nolep_30",2200,0,"M_{J} without leptons",2,1));
  mj_lep.push_back(var_class("mj_siglep_30",2200,0,"M_{J} with 20 GeV leptons",28,1));
  mj_lep.push_back(var_class("mj_30", 2200,0,"M_{J} with 30 GeV leptons",4,1));

  ///////////////////// ROCs to plot /////////////////////
  vector<TString> vs_sam, vs_vars;
  vector<vector<sample_class>*> v_sam; 
  v_sam.push_back(&qcd_t1t); vs_sam.push_back("qcd_t1t");
  v_sam.push_back(&qcd_t1tc); vs_sam.push_back("qcd_t1tc");
  v_sam.push_back(&qcd_t1q); vs_sam.push_back("qcd_t1q");
  v_sam.push_back(&qcd_t1qc); vs_sam.push_back("qcd_t1qc");
  v_sam.push_back(&tt_t1t); vs_sam.push_back("tt_t1t");
  v_sam.push_back(&tt_t1tc); vs_sam.push_back("tt_t1tc");
  v_sam.push_back(&tt_t1q); vs_sam.push_back("tt_t1q");
  v_sam.push_back(&tt_t1qc); vs_sam.push_back("tt_t1qc");
  v_sam.push_back(&zjets_t1t); vs_sam.push_back("zjets_t1t");
  v_sam.push_back(&zjets_t1tc); vs_sam.push_back("zjets_t1tc");
  v_sam.push_back(&zjets_t1q); vs_sam.push_back("zjets_t1q");
  v_sam.push_back(&zjets_t1qc); vs_sam.push_back("zjets_t1qc");

  vector<vector<var_class>*> v_vars;
  // v_vars.push_back(&mj_general); vs_vars.push_back("general");
  // v_vars.push_back(&mj_sizes); vs_vars.push_back("size");
  // v_vars.push_back(&mj_cands); vs_vars.push_back("cands");
  // v_vars.push_back(&mj_pt); vs_vars.push_back("pt");
  v_vars.push_back(&mj_eta); vs_vars.push_back("eta"); 
  v_vars.push_back(&mj_ptfat); vs_vars.push_back("ptfat"); 

  vector<TString> v_cuts;
  v_cuts.push_back("ht>500&&met>200&&njets30>=4&&nvmus10==0&&nvels10==0");

  for(unsigned ivar(0); ivar<v_vars.size(); ivar++){
    for(unsigned icut(0); icut<v_cuts.size(); icut++){
      for(unsigned isam(0); isam<v_sam.size(); isam++){
  	DrawROC(*(v_sam[isam]), *(v_vars[ivar]), v_cuts[icut], "mj_"+vs_sam[isam]+"_"+vs_vars[ivar]);
      } // Loop over samples
    } // Loop over cuts
  } // Loop over variables

  v_sam.clear(); vs_sam.clear();
  v_sam.push_back(&tt_t1t); vs_sam.push_back("tt_t1t");
  v_sam.push_back(&tt_t1tc); vs_sam.push_back("tt_t1tc");

  v_vars.push_back(&mj_lep); vs_vars.push_back("lep"); 

  v_cuts.clear();
  v_cuts.push_back("ht>500&&met>200&&njets>=4&&(nmus+nels)==1");
  v_cuts.push_back("ht>500&&met>200&&njets>=4&&(nmus+nels)==2");

  for(unsigned ivar(0); ivar<v_vars.size(); ivar++){
    for(unsigned icut(0); icut<v_cuts.size(); icut++){
      for(unsigned isam(0); isam<v_sam.size(); isam++){
  	DrawROC(*(v_sam[isam]), *(v_vars[ivar]), v_cuts[icut], "mj_"+vs_sam[isam]+"_"+vs_vars[ivar]);
      } // Loop over samples
    } // Loop over cuts
  } // Loop over variables

  //DrawROC(tt_t1t, mj_lep, "ht>500&&met>200&&njets30>=4&&(nmus+nels)==1", "mj_tt_lep");
  //DrawROC(tt_t1t, mj_lep, "ht>500&&met>200&&njets30>=4&&(nmus+nels)==2", "mj_tt_lep");
}



void DrawROC(vector<sample_class> samples, vector<var_class> vars, TString cuts, TString tag){
  TCanvas can;
  const int nbins(1000);
  vector<vector<TH1D> > histos;
  TString hname, totcut;
  TChain *chain[2];

  for(unsigned sam(0); sam < samples.size(); sam++){
    // Loading chains
    for(unsigned isam(0); isam < samples[sam].files.size(); isam++){
      chain[sam] = new TChain("tree");
      int nfiles = chain[sam]->Add(samples[sam].files[isam]);
      if(nfiles==0) cout<<samples[sam].files[isam]<<" not found"<<endl;
    }
    histos.push_back(vector<TH1D>());

    // Projecting variables
    for(unsigned var(0); var<vars.size(); var++){
      float minh(vars[var].minx), maxh(vars[var].maxx);
      if(minh > maxh){
	minh = maxh;
	maxh = vars[var].minx;
      }
      hname = "histo"; hname += sam; hname += var;
      totcut = "weight*("+cuts+"&&"+samples[sam].cut+")";
      histos[sam].push_back(TH1D(hname,"",nbins,minh,maxh));
      chain[sam]->Project(hname, vars[var].varname, totcut);
    } // Loop over variables
  } // Loop over samples

  TString title(cuts);
  if(title=="1") title = "";
  title.ReplaceAll("&&1",""); title.ReplaceAll("nvmus10==0&&nvels10==0", "0 leptons");  
  title.ReplaceAll("(nmus+nels)", "n_{lep}");  
  title.ReplaceAll("els_pt","p^{e}_{T}");title.ReplaceAll("mus_pt","p^{#mu}_{T}");
  title.ReplaceAll("mj_30", "M_{J}");
  title.ReplaceAll(">=", " #geq "); 
  title.ReplaceAll(">", " > "); title.ReplaceAll("<", " < "); title.ReplaceAll("&&", ", "); 
  title.ReplaceAll("met", "MET");  title.ReplaceAll("ht30", "H_{T}^{30}"); title.ReplaceAll("ht", "H_{T}"); 
  title.ReplaceAll("mt", "m_{T}"); 
  title.ReplaceAll("nleps==1", "1 lepton");  title.ReplaceAll("ntrupv","n^{true}_{PV}");
  title.ReplaceAll("njets30","n_{jets}^{30}"); title.ReplaceAll("nbm30","n_{b}^{30}");
  title.ReplaceAll("njets","n_{jets}"); title.ReplaceAll("nbm","n_{b}");
  title.ReplaceAll("mindphin_metje","min#Delta#phi_{N}");
  title.ReplaceAll("nbl","n_{b,l}");
  TH1D base_histo("base",title,1,0.01,1.0);
  base_histo.SetXTitle(samples[0].label+" efficiency");
  base_histo.SetYTitle(samples[1].label+" efficiency");
  base_histo.SetMinimum(0.0);
  base_histo.SetMaximum(1.0);
  base_histo.SetDirectory(0);
  base_histo.Draw();

  // Legend
  double legX = 0.13, legY = 0.88, legSingle = 0.064;
  double legW = 0.12, legH = legSingle*vars.size();
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(0.055); leg.SetFillColor(0); leg.SetFillStyle(0); leg.SetBorderSize(0);
  leg.SetTextFont(132);

  // Making individual graphs
  TGraph graphs[100]; // Had to make it an array because a vector<TGraph> kept crashing
  for(unsigned var(0); var<vars.size(); var++){
    graphs[var] = MakeROC(histos[0][var], histos[1][var], vars[var].minx < vars[var].maxx, vars[var].cuts);
    graphs[var].SetLineColor(vars[var].color);
    graphs[var].SetLineStyle(vars[var].style);
    graphs[var].SetLineWidth(3);
    leg.AddEntry(&(graphs[var]), vars[var].title, "l");
    graphs[var].Draw("lsame");
  } // Loop over variables
  leg.Draw();

  cuts.ReplaceAll(".",""); 
  cuts.ReplaceAll("(",""); cuts.ReplaceAll("$","");  cuts.ReplaceAll(")",""); 
  cuts.ReplaceAll("[",""); cuts.ReplaceAll("]",""); 
  cuts.ReplaceAll("/","_"); cuts.ReplaceAll("*",""); cuts.ReplaceAll("&&","_");
  cuts.ReplaceAll(">=","ge"); cuts.ReplaceAll("<=","se"); 
  cuts.ReplaceAll(">","g"); cuts.ReplaceAll("<","s"); cuts.ReplaceAll("=","");
  cuts.ReplaceAll("+",""); 
  TString pname("plots/roc_"+tag+"_"+cuts+".pdf");  
  can.Print(pname);
  can.SetLogx(1);
  can.SetLogy(1);
  pname.ReplaceAll("roc_","log_roc_");
  base_histo.SetMinimum(1e-4);
  can.Print(pname);

  for(unsigned sam(0); sam < samples.size(); sam++)
    chain[sam]->Delete();
}

TGraph MakeROC(TH1D &good, TH1D &bad, const bool less_is_better, vector<marker_class> cuts){
  const int nbins = good.GetNbinsX();
  if(bad.GetNbinsX() != nbins) throw logic_error("Mismatched number of bins");

  TMarker marker;  

  TGraph graph(0);
  const double good_tot = good.Integral(0, nbins+1);
  const double bad_tot = bad.Integral(0, nbins+1);
  int inibin(0), endbin(nbins+1), dbin(1); unsigned icut(0);
  if(less_is_better){
    inibin = nbins+1;
    endbin = 0;
    dbin = -1;
  }
  for(int bin = inibin; bin*dbin<=endbin*dbin; bin+=dbin){
    const double good_pass = good.Integral(min(endbin,bin), max(endbin,bin));
    const double bad_pass = bad.Integral(min(endbin,bin), max(endbin,bin));
    const double x = good_pass/good_tot;
    const double y = bad_pass/bad_tot;
    graph.SetPoint(graph.GetN(), x, y);
 
    // Plotting the stars
    if(icut<cuts.size()){
      float edge(good.GetXaxis()->GetBinUpEdge(bin));
      if((edge>=cuts[icut].cut&&!less_is_better) || (edge<=cuts[icut].cut&&less_is_better)){
	marker.SetMarkerStyle(cuts[icut].style);marker.SetMarkerColor(cuts[icut].color);
	marker.SetMarkerSize(cuts[icut].size); 
	marker.DrawMarker(x,y);
	icut++;
      }
    }
  }
  TString name(good.GetName());
  name += "graph";
  graph.SetName(name);
  graph.SetTitle(name);

  return graph;
}

var_class::var_class(TString ivarname, float iminx, float imaxx, TString ititle, int icolor, 
	    int istyle, vector<marker_class> icuts){
  varname = ivarname; minx = iminx; maxx = imaxx; title = ititle;
  cuts = icuts; 
  color = icolor; style = istyle;
}

sample_class::sample_class(TString ilabel, vector<TString> ifiles, TString icut){
  files = ifiles; label = ilabel; cut = icut;
}

marker_class::marker_class(float icut, float isize, int icolor, int istyle){
  cut=icut; size=isize; color=icolor; style=istyle;
}

