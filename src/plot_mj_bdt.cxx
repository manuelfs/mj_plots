// plot_mj_bdt: Macro that plots ROC curves comparing the single variable performance
//              to various BDTs

#include <stdexcept>
#include <iostream>

#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TDirectory.h"
#include "TMarker.h"
#include "TStyle.h"

#include "styles.hpp"
#include "utilities.hpp"
#include "plot_mj_bdt.hpp"
#include "utilities_macros.hpp"

using namespace std;

int main(){
  styles style("Paper"); style.setDefaultStyle();
  gStyle->SetPadTickX(1); // Tickmarks on the top
  gStyle->SetPadTickY(1); // Tickmarks on the right
  TString folder="root/";

  ////////////////////////// SAMPLES for the axes /////////////////////////
  vector<TString> v_t1t_tt_htmj_30;
  v_t1t_tt_htmj_30.push_back(folder+"bdt_T1tttt1500_TTJet_htmj_30.root");
  vector<sample_class> t1t_tt_htmj_30; 
  t1t_tt_htmj_30.push_back(sample_class("T1tttt(1500,100)", v_t1t_tt_htmj_30));
  t1t_tt_htmj_30.push_back(sample_class("t#bar{t}", v_t1t_tt_htmj_30));

  vector<TString> v_t1t_tt_htnjets30;
  v_t1t_tt_htnjets30.push_back(folder+"bdt_T1tttt1500_TTJet_htnjets30.root");
  vector<sample_class> t1t_tt_htnjets30; 
  t1t_tt_htnjets30.push_back(sample_class("T1tttt(1500,100)", v_t1t_tt_htnjets30));
  t1t_tt_htnjets30.push_back(sample_class("t#bar{t}", v_t1t_tt_htnjets30));

  vector<TString> v_t1t_tt_njets30mj_30;
  v_t1t_tt_njets30mj_30.push_back(folder+"bdt_T1tttt1500_TTJet_njets30mj_30.root");
  vector<sample_class> t1t_tt_njets30mj_30; 
  t1t_tt_njets30mj_30.push_back(sample_class("T1tttt(1500,100)", v_t1t_tt_njets30mj_30));
  t1t_tt_njets30mj_30.push_back(sample_class("t#bar{t}", v_t1t_tt_njets30mj_30));

  // T1tttt(1200,800)
  vector<TString> v_t1tc_tt_htmj_30;
  v_t1tc_tt_htmj_30.push_back(folder+"bdt_T1tttt1200_TTJet_htmj_30.root");
  vector<sample_class> t1tc_tt_htmj_30; 
  t1tc_tt_htmj_30.push_back(sample_class("T1tttt(1200,800)", v_t1tc_tt_htmj_30));
  t1tc_tt_htmj_30.push_back(sample_class("t#bar{t}", v_t1tc_tt_htmj_30));

  vector<TString> v_t1tc_tt_htnjets30;
  v_t1tc_tt_htnjets30.push_back(folder+"bdt_T1tttt1500_TTJet_htnjets30.root");
  vector<sample_class> t1tc_tt_htnjets30; 
  t1tc_tt_htnjets30.push_back(sample_class("T1tttt(1500,100)", v_t1tc_tt_htnjets30));
  t1tc_tt_htnjets30.push_back(sample_class("t#bar{t}", v_t1tc_tt_htnjets30));

  vector<TString> v_t1tc_tt_njets30mj_30;
  v_t1tc_tt_njets30mj_30.push_back(folder+"bdt_T1tttt1500_TTJet_njets30mj_30.root");
  vector<sample_class> t1tc_tt_njets30mj_30; 
  t1tc_tt_njets30mj_30.push_back(sample_class("T1tttt(1500,100)", v_t1tc_tt_njets30mj_30));
  t1tc_tt_njets30mj_30.push_back(sample_class("t#bar{t}", v_t1tc_tt_njets30mj_30));

  // T1qqqq(1400,100)
  vector<TString> v_t1q_tt_htmj_30;
  v_t1q_tt_htmj_30.push_back(folder+"bdt_T1qqqq1400_TTJet_htmj_30.root");
  vector<sample_class> t1q_tt_htmj_30; 
  t1q_tt_htmj_30.push_back(sample_class("T1qqqq(1400,100)", v_t1q_tt_htmj_30));
  t1q_tt_htmj_30.push_back(sample_class("t#bar{t}", v_t1q_tt_htmj_30));

  vector<TString> v_t1q_tt_htnjets30;
  v_t1q_tt_htnjets30.push_back(folder+"bdt_T1qqqq1400_TTJet_htnjets30.root");
  vector<sample_class> t1q_tt_htnjets30; 
  t1q_tt_htnjets30.push_back(sample_class("T1qqqq(1400,100)", v_t1q_tt_htnjets30));
  t1q_tt_htnjets30.push_back(sample_class("t#bar{t}", v_t1q_tt_htnjets30));

  vector<TString> v_t1q_tt_njets30mj_30;
  v_t1q_tt_njets30mj_30.push_back(folder+"bdt_T1qqqq1400_TTJet_njets30mj_30.root");
  vector<sample_class> t1q_tt_njets30mj_30; 
  t1q_tt_njets30mj_30.push_back(sample_class("T1qqqq(1400,100)", v_t1q_tt_njets30mj_30));
  t1q_tt_njets30mj_30.push_back(sample_class("t#bar{t}", v_t1q_tt_njets30mj_30));


  ///////////////////// Markers for each ROC /////////////////////
  int ht_col(2);
  int mj_style(8); float mj_size(2.5);
  vector<marker_class> mj_points, ht_points, nj_points;
  mj_points.push_back(marker_class(200,  mj_size, 4, mj_style));
  mj_points.push_back(marker_class(400,  mj_size, 4, mj_style));
  mj_points.push_back(marker_class(600,  4, 4, 29));
  mj_points.push_back(marker_class(800,  mj_size, 4, mj_style));
  mj_points.push_back(marker_class(1000,  mj_size, 4, mj_style));
  ht_points.push_back(marker_class(1000,  mj_size, ht_col, mj_style));
  ht_points.push_back(marker_class(1500,  4, ht_col, 29));
  ht_points.push_back(marker_class(2000,  mj_size, ht_col, mj_style));
  ht_points.push_back(marker_class(2500,  mj_size, ht_col, mj_style));
  nj_points.push_back(marker_class(6,  mj_size, 28, mj_style));
  nj_points.push_back(marker_class(8,  4, 28, 29));
  nj_points.push_back(marker_class(10,  mj_size, 28, mj_style));
  nj_points.push_back(marker_class(12,  mj_size, 28, mj_style));


  ///////////////////// ROCs to plot /////////////////////
  vector<var_class> mj_t1t_tt;
  mj_t1t_tt.push_back(var_class(t1t_tt_htmj_30, "ht",4000,0,"H_{T}",ht_col,1,ht_points));
  mj_t1t_tt.push_back(var_class(t1t_tt_htmj_30, "njets30",15,0,"n_{jets}^{30}",28,1,nj_points));
  mj_t1t_tt.push_back(var_class(t1t_tt_htmj_30, "mj_30",2000,0,"M_{J}",4,1,mj_points));
  // mj_t1t_tt.push_back(var_class(t1t_tt_htmj_30, "BDT",0.22, -0.5,"BDT[H_{T}, M_{J}]",kMagenta+1,7));
  // mj_t1t_tt.push_back(var_class(t1t_tt_htnjets30, "BDT",0.3, -0.5,"BDT[H_{T}, n_{jets}^{30}]",1,2));
  // mj_t1t_tt.push_back(var_class(t1t_tt_njets30mj_30, "BDT",0.32, -0.5,"BDT[M_{J}, n_{jets}^{30}]",kGreen+2,2));

  vector<var_class> mj_t1tc_tt;
  mj_t1tc_tt.push_back(var_class(t1tc_tt_htmj_30, "ht",4000,0,"H_{T}",ht_col,1,ht_points));
  mj_t1tc_tt.push_back(var_class(t1tc_tt_htmj_30, "njets30",15,0,"n_{jets}^{30}",28,1,nj_points));
  mj_t1tc_tt.push_back(var_class(t1tc_tt_htmj_30, "mj_30",2000,0,"M_{J}",4,1,mj_points));
  mj_t1tc_tt.push_back(var_class(t1tc_tt_htmj_30, "BDT",0.22, -0.5,"BDT[H_{T}, M_{J}]",kMagenta+1,7));
  mj_t1tc_tt.push_back(var_class(t1tc_tt_htnjets30, "BDT",0.3, -0.5,"BDT[H_{T}, n_{jets}^{30}]",1,2));
  mj_t1tc_tt.push_back(var_class(t1tc_tt_njets30mj_30, "BDT",0.32, -0.5,"BDT[M_{J}, n_{jets}^{30}]",kGreen+2,2));

  vector<var_class> mj_t1q_tt;
  mj_t1q_tt.push_back(var_class(t1q_tt_htmj_30, "ht",4000,0,"H_{T}",ht_col,1,ht_points));
  mj_t1q_tt.push_back(var_class(t1q_tt_htmj_30, "njets30",15,0,"n_{jets}^{30}",28,1,nj_points));
  mj_t1q_tt.push_back(var_class(t1q_tt_htmj_30, "mj_30",2000,0,"M_{J}",4,1,mj_points));
  mj_t1q_tt.push_back(var_class(t1q_tt_htmj_30, "BDT",0.22, -0.5,"BDT[H_{T}, M_{J}]",kMagenta+1,7));
  mj_t1q_tt.push_back(var_class(t1q_tt_htnjets30, "BDT",0.3, -0.5,"BDT[H_{T}, n_{jets}^{30}]",1,2));
  mj_t1q_tt.push_back(var_class(t1q_tt_njets30mj_30, "BDT",0.32, -0.5,"BDT[M_{J}, n_{jets}^{30}]",kGreen+2,2));

  vector<TString> vs_sam;
  vector<vector<var_class> > v_mj;
  v_mj.push_back(mj_t1t_tt);   vs_sam.push_back("t1t_tt");
  v_mj.push_back(mj_t1tc_tt);  vs_sam.push_back("t1tc_tt");
  v_mj.push_back(mj_t1q_tt);   vs_sam.push_back("t1q_tt");

  vector<TString> v_cuts;
  v_cuts.push_back("ht>500&&met>200&&njets30>=4&&nvmus10==0&&nvels10==0");
  //v_cuts.push_back("ht>500&&met>200&&met<=400&&njets30>=4&&nvmus10==0&&nvels10==0");
  //v_cuts.push_back("ht>500&&met>400&&njets30>=4&&nvmus10==0&&nvels10==0");

  for(unsigned icut(0); icut<v_cuts.size(); icut++){
    for(unsigned isam(0); isam<vs_sam.size(); isam++){
      DrawROC(v_mj[isam], v_cuts[icut], "mj_nobdt_"+vs_sam[isam]);
    } // Loop over samples
  } // Loop over cuts

}



void DrawROC(vector<var_class> vars, TString cuts, TString tag){
  TCanvas can;
  const int nbins(1000);
  vector<vector<TH1D> > histos;
  TString hname, totcut;
  TChain *chain[2];

  for(unsigned var(0); var<vars.size(); var++){
    vector<sample_class> samples = vars[var].samples;
    for(unsigned sam(0); sam < samples.size(); sam++){
      // Loading chains
      for(unsigned isam(0); isam < samples[sam].files.size(); isam++){
	chain[sam] = new TChain("TestTree");
	int nfiles = chain[sam]->Add(samples[sam].files[isam]);
	if(nfiles==0) cout<<samples[sam].files[isam]<<" not found"<<endl;
      }
      histos.push_back(vector<TH1D>());

      // Projecting variables
      float minh(vars[var].minx), maxh(vars[var].maxx);
      if(minh > maxh){
	minh = maxh;
	maxh = vars[var].minx;
      }
      hname = "histo"; hname += sam; hname += var;
      totcut = "weight*("+cuts+"&&"+samples[sam].cut+"&&classID==";
      totcut += sam; totcut += ")";
      histos[sam].push_back(TH1D(hname,"",nbins,minh,maxh));
      chain[sam]->Project(hname, vars[var].varname, totcut);
      chain[sam]->Delete();
    } // Loop over variables
  } // Loop over samples

  TString title(cuts);
  if(title=="1") title = "";
  title.ReplaceAll("&&1",""); title.ReplaceAll("nvmus10==0&&nvels10==0", "0 leptons");  
  title.ReplaceAll("(nmus+nels)", "n_{lep}");  
  title.ReplaceAll("els_pt","p^{e}_{T}");title.ReplaceAll("mus_pt","p^{#mu}_{T}");
  title.ReplaceAll("mj_30", "M_{J}"); title.ReplaceAll("<=", " #leq "); 
  title.ReplaceAll(">=", " #geq "); title.ReplaceAll("==", " = "); 
  title.ReplaceAll(">", " > "); title.ReplaceAll("<", " < "); title.ReplaceAll("&&", ", "); 
  title.ReplaceAll("met", "MET");  title.ReplaceAll("ht30", "H_{T}^{30}"); title.ReplaceAll("ht", "H_{T}"); 
  title.ReplaceAll("mt", "m_{T}"); 
  title.ReplaceAll("nleps==1", "1 lepton");  title.ReplaceAll("ntrupv","n^{true}_{PV}");
  title.ReplaceAll("njets30","n_{jets}^{30}"); title.ReplaceAll("nbm30","n_{b}^{30}");
  title.ReplaceAll("njets","n_{jets}"); title.ReplaceAll("nbm","n_{b}");
  title.ReplaceAll("mindphin_metje","min#Delta#phi_{N}");
  title.ReplaceAll("nbl","n_{b,l}");
  TH1D base_histo("base",title,1,0.01,1.0);
  base_histo.SetXTitle(vars[0].samples[0].label+" efficiency");
  base_histo.SetYTitle(vars[0].samples[1].label+" efficiency");
  base_histo.SetMinimum(0.0);
  base_histo.SetMaximum(1.0);
  base_histo.SetDirectory(0);
  base_histo.Draw();

  // Legend
  double legX = 0.14, legY = 0.88, legSingle = 0.081;
  double legW = 0.2, legH = legSingle*vars.size();
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(0.062); leg.SetFillColor(0); leg.SetFillStyle(0); leg.SetBorderSize(0);
  leg.SetTextFont(132);

  // Making individual graphs
  TGraph graphs[100]; // Had to make it an array because a vector<TGraph> kept crashing
  for(unsigned var(0); var<vars.size(); var++){
    graphs[var] = MakeROC(histos[0][var], histos[1][var], vars[var].minx < vars[var].maxx, vars[var].cuts);
    graphs[var].SetLineColor(vars[var].color);
    graphs[var].SetLineStyle(vars[var].style);
    if(vars[var].style==1) graphs[var].SetLineWidth(5);
    else  graphs[var].SetLineWidth(8);
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
	if(x>0.01 && y>0.0001) marker.DrawMarker(x,y);
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

var_class::var_class(vector<sample_class> isamples, TString ivarname, float iminx, float imaxx, TString ititle, 
		     int icolor, int istyle, vector<marker_class> icuts){
  varname = ivarname; minx = iminx; maxx = imaxx; title = ititle;
  cuts = icuts; 
  color = icolor; style = istyle;
  samples = isamples;
}

sample_class::sample_class(TString ilabel, vector<TString> ifiles, TString icut){
  files = ifiles; label = ilabel; cut = icut;
}

marker_class::marker_class(float icut, float isize, int icolor, int istyle){
  cut=icut; size=isize; color=icolor; style=istyle;
}

