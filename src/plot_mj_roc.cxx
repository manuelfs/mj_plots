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

  TString folder="/cms2r0/manuelf/root/small/archive/15-01-08/";

  // NTUPLES
  vector<TString> v_t1;
  v_t1.push_back(folder+"*T1tttt*1500_*PU20*");

  vector<TString> v_tt;
  v_tt.push_back(folder+"/skim/small_TTJet_ht30_500_met200.root");

  vector<TString> v_qcd;
  v_qcd.push_back(folder+"/skim/small_qcd_ht30_500_met200.root");

  ////////////////////////// SAMPLES for the axes /////////////////////////
  vector<sample_class> tt_t1; 
  tt_t1.push_back(sample_class("T1tttt(1500,100)", v_t1));
  tt_t1.push_back(sample_class("t#bar{t}", v_tt));

  vector<sample_class> qcd_t1; 
  qcd_t1.push_back(sample_class("T1tttt(1500,100)", v_t1));
  qcd_t1.push_back(sample_class("QCD", v_qcd));



  ///////////////////// VARIABLES for each ROC /////////////////////
  vector<var_class> mj_vars;
  mj_vars.push_back(var_class("ht30",4000,0,"H_{T}^{30}",4,1));
  mj_vars.push_back(var_class("met",1500,0,"MET",2,1));
  mj_vars.push_back(var_class("mj_30",2000,0,"M_{J}",8,1));

  vector<var_class> mj_sizes;
  mj_sizes.push_back(var_class("ht30",4000,0,"H_{T}^{30}",2,1));
  mj_sizes.push_back(var_class("mj_r08",2200,0,"M_{J} R=0.8",8,1));
  mj_sizes.push_back(var_class("mj_r10",2200,0,"M_{J} R=1.0",28,1));
  mj_sizes.push_back(var_class("mj_30", 2200,0,"M_{J} R=1.2",4,1));
  mj_sizes.push_back(var_class("mj_r14",2200,0,"M_{J} R=1.4",kMagenta+2,1));

  vector<var_class> mj_cands;
  mj_cands.push_back(var_class("ht30",4000,0,"H_{T}^{30}",2,1));
  mj_cands.push_back(var_class("mj_cands",2200,0,"M_{J} pfcands",8,1));
  mj_cands.push_back(var_class("mj_cands_trim",2200,0,"M_{J} pfcands trimmed",28,1));
  mj_cands.push_back(var_class("mj_30", 2200,0,"M_{J} 30 GeV jets",4,1));

  vector<var_class> mj_pt;
  mj_pt.push_back(var_class("ht30",4000,0,"H_{T}^{30}",2,1));
  mj_pt.push_back(var_class("mj_10",2200,0,"M_{J} 10 GeV jets",8,1));
  mj_pt.push_back(var_class("mj_20",2200,0,"M_{J} 20 GeV jets",28,1));
  mj_pt.push_back(var_class("mj_30", 2200,0,"M_{J} 30 GeV jets",4,1));
  mj_pt.push_back(var_class("mj_40",2200,0,"M_{J} 40 GeV jets",kMagenta+2,1));

  ///////////////////// ROCs to plot /////////////////////
  vector<TString> vs_sam, vs_vars;
  vector<vector<sample_class>*> v_sam; 
  v_sam.push_back(&tt_t1); vs_sam.push_back("tt");
  v_sam.push_back(&qcd_t1); vs_sam.push_back("qcd");

  vector<vector<var_class>*> v_vars;
  v_vars.push_back(&mj_vars); vs_vars.push_back("general");
  v_vars.push_back(&mj_sizes); vs_vars.push_back("size");
  v_vars.push_back(&mj_cands); vs_vars.push_back("cands");
  v_vars.push_back(&mj_pt); vs_vars.push_back("pt");

  TString cuts("ht30>500&&met>200&&njets30>=4");
  for(unsigned ivar(0); ivar<v_vars.size(); ivar++){
    for(unsigned isam(0); isam<v_sam.size(); isam++){
      DrawROC(*(v_sam[isam]), *(v_vars[ivar]), cuts, "mj_"+vs_sam[isam]+"_"+vs_vars[ivar]);
    }
  }
  // DrawROC(tt_t1, mj_vars, "ht30>500&&met>200&&njets30>=4", "mj_tt");
  // DrawROC(qcd_t1, mj_vars, "ht30>500&&met>200&&njets30>=4", "mj_qcd");
  // DrawROC(tt_t1, mj_sizes, "ht30>500&&met>200&&njets30>=4", "mj_tt_sizes");
  // DrawROC(tt_t1, mj_cands, "ht30>500&&met>200&&njets30>=4", "mj_tt_cands");
  // DrawROC(qcd_t1, mj_sizes, "ht30>500&&met>200&&njets30>=4", "mj_qcd_sizes");
  // DrawROC(qcd_t1, mj_cands, "ht30>500&&met>200&&njets30>=4", "mj_qcd_cands");
  // DrawROC(tt_t1, mj_pt, "ht30>500&&met>200&&njets30>=4", "mj_tt_pt");

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
  TString pname("eps/roc_"+tag+"_"+cuts+".eps");  
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

