// plot_distribution: Macro that plots variables both lumi weighted and normalized to the same area.

#include <iostream>
#include <vector>

#include "TChain.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"
#include "TColor.h"
#include "TF1.h"

#include "styles.hpp"
#include "utilities.hpp"
#include "small_tree.hpp"

using namespace std;
using std::cout;
using std::endl;

class hfeats {
public:
  hfeats(TString ivarname, int inbins, float iminx, float imaxx, vector<int> isamples,
	 TString ititle="", TString icuts="1", float icut=-1){
    varname = ivarname; nbins = inbins; minx = iminx; maxx = imaxx; title = ititle;
    cuts = icuts; cut = icut; samples = isamples;
    tag = ivarname+"_"+cuts; tag.ReplaceAll("_1",""); tag.ReplaceAll(".",""); tag.ReplaceAll(" ","_"); 
    tag.ReplaceAll("(",""); tag.ReplaceAll("$","");  tag.ReplaceAll(")",""); 
    tag.ReplaceAll("[",""); tag.ReplaceAll("]",""); tag.ReplaceAll("||","_");
    tag.ReplaceAll("/","_"); tag.ReplaceAll("*",""); tag.ReplaceAll("&&","_");
    tag.ReplaceAll(">",""); tag.ReplaceAll("<",""); tag.ReplaceAll("=","");
    tag.ReplaceAll("+",""); tag.ReplaceAll("{",""); tag.ReplaceAll("}",""); 
    tag.ReplaceAll("^",""); tag.ReplaceAll(",",""); 
    unit = "";
    if(title.Contains("GeV)")) unit = "GeV";
    if(title.Contains("phi")) unit = "rad";
  }
  TString title, varname, tag, cuts, unit;
  int nbins;
  float minx, maxx, cut;
  vector<int> samples;
};

class sfeats {
public:
  sfeats(vector<TString> ifile, TString ilabel, int icolor, int istyle=1, TString icut="1"){
    file = ifile; label = ilabel; cut = icut;
    color = icolor; style = istyle;
    isSig = ifile[0].Contains("T1tttt");// && ifile.Contains("1200");
    factor = "1";
    if(ifile[0].Contains("TTW")) factor = "0.36";
  }
  vector<TString> file;
  TString label, cut, factor;
  int color, style;
  bool isSig;
};

int main(){ 
  styles style("Paper"); style.PadLeftMargin = 0.13; style.yTitleOffset = 1.07;
  style.setDefaultStyle(); 
  vector<hfeats> vars;
  TCanvas can;

  TColor ucsb_blue(1000, 1/255.,57/255.,166/255.);
  TColor ucsb_gold(1001, 255/255.,200/255.,47/255);
  TColor penn_red(1002, 149/255.,0/255.,26/255.);
  TColor uf_orange(1003, 255/255.,74/255.,0/255.);
  TColor uo_green(1004, 0/255.,79/255.,39/255.);
  TColor tcu_purple(1005, 52/255.,42/255.,123/255.);
  TColor tar_heel_blue(1006, 86/255.,160/255.,211/255.);
  TColor sig_teal(1007, 96/255.,159/255.,128/255.);
  TColor sig_gold(1008, 215/255.,162/255.,50/255.);
  TColor seal_brown(1010, 89/255.,38/255.,11/255.);

  //TString folder="archive/15-03-03/skim/";
  TString folder="archive/15-01-30/skim/";
  vector<TString> s_tt;
  s_tt.push_back(folder+"*_TTJet*");
  vector<TString> s_wjets;
  s_wjets.push_back(folder+"*WJets*");
  vector<TString> s_zjets;
  s_zjets.push_back(folder+"*_ZJets*");
  vector<TString> s_single;
  s_single.push_back(folder+"*_T*channel*");
  vector<TString> s_ttv;
  s_ttv.push_back(folder+"*TTW*");
  s_ttv.push_back(folder+"*TTZ*");
  vector<TString> s_qcd;
  s_qcd.push_back(folder+"*QCD_Pt*");
  vector<TString> s_other;
  s_other.push_back(folder+"*DY*");
  s_other.push_back(folder+"*WH_HToBB*");
  vector<TString> s_t1t;
  s_t1t.push_back(folder+"*T1tttt*1500_*PU20*");
  vector<TString> s_t1tc;
  s_t1tc.push_back(folder+"*T1tttt*1200_*PU20*");
  vector<TString> s_t1q;
  s_t1q.push_back(folder+"*T1qqqq*1400_*PU20*");
  vector<TString> s_t1qc;
  s_t1qc.push_back(folder+"*T1qqqq*1000_*PU20*");

  // Reading ntuples
  vector<TChain *> chain;
  vector<sfeats> Samples; 
  // MJ note's colors
  Samples.push_back(sfeats(s_t1t, "T1tttt(1500,100)", 4));
  Samples.push_back(sfeats(s_t1q, "T1qqqq(1400,100)", 8));
  Samples.push_back(sfeats(s_zjets, "Z#rightarrow#nu#nu", kMagenta+2));
  Samples.push_back(sfeats(s_qcd, "QCD", 28));
  Samples.push_back(sfeats(s_tt, "t#bar{t}", 2));
  //Samples.push_back(sfeats(s_t1qc, "T1qqqq(1000,800)", 8,2));
  //Samples.push_back(sfeats(s_t1tc, "T1tttt(1200,800)", 4,2));

  // // Jack's colors
  // Samples.push_back(sfeats(s_zjets, "Z#rightarrow#nu#nu", 1002));
  // Samples.push_back(sfeats(s_qcd, "QCD", 1001));
  // Samples.push_back(sfeats(s_tt, "t#bar{t}", 1000,1));
  // Samples.push_back(sfeats(s_t1q, "T1qqqq(1400,100)", 4));
  // Samples.push_back(sfeats(s_t1t, "T1tttt(1500,100)", 2));
  // Samples.push_back(sfeats(s_t1qc, "T1qqqq(1000,800)", 4,2));
  // Samples.push_back(sfeats(s_t1tc, "T1tttt(1200,800)", 2,2));

  // Samples.push_back(sfeats(s_ttv, "ttV", 1002));
  // Samples.push_back(sfeats(s_single, "Single top", 1005));
  // Samples.push_back(sfeats(s_wjets, "W + jets", 1004));

  for(unsigned sam(0); sam < Samples.size(); sam++){
    chain.push_back(new TChain("tree"));
    for(unsigned insam(0); insam < Samples[sam].file.size(); insam++)
      chain[sam]->Add(Samples[sam].file[insam]);
  }

  vector<int> mj_sam;
  mj_sam.push_back(0);
  mj_sam.push_back(1);
  mj_sam.push_back(2);
  mj_sam.push_back(3);
  mj_sam.push_back(4);

  vars.push_back(hfeats("Average M_{J}^{10} (GeV)",7,5,40, mj_sam, "True n_{PV}",
   			"ht>500&&met>200&&njets>=4&&(nmus+nels)==0"));
  vars.push_back(hfeats("Average M_{J}^{30} (GeV)",7,5,40, mj_sam, "True n_{PV}",
   			"ht>500&&met>200&&njets>=4&&(nmus+nels)==0"));
  vars.push_back(hfeats("Average M_{J}^{40} (GeV)",7,5,40, mj_sam, "True n_{PV}",
   			"ht>500&&met>200&&njets>=4&&(nmus+nels)==0"));
  vars.push_back(hfeats("Average M_{J}^{pfcand} (GeV)",7,5,40, mj_sam, "True n_{PV}",
			"ht>500&&met>200&&njets>=4&&(nmus+nels)==0"));
  vars.push_back(hfeats("Average M_{J}^{pfcand,trim} (GeV)",7,5,40, mj_sam, "True n_{PV}",
			"ht>500&&met>200&&njets>=4&&(nmus+nels)==0"));

  double legX = 0.15, legY = 0.89, legSingle = 0.067;
  double legW = 0.22, legH = legSingle*Samples.size();
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(0.057); leg.SetFillColor(0); leg.SetFillStyle(0); leg.SetBorderSize(0);
  leg.SetTextFont(132);

  vector<vector<TH1D*> > varhisto[4];
  TString hname, pname, leghisto, title;
  TF1 linfit("linfit","[0]+[1]*x",vars[0].minx, vars[0].maxx);
  linfit.SetParameter(0, 200);  linfit.SetParameter(1, 0);

  for(unsigned sam(0); sam < Samples.size(); sam++){

    // Setting up histogramms
    for(unsigned his(0); his < 4; his++){    
      varhisto[his].push_back(vector<TH1D*>());
      for(unsigned var(0); var<vars.size(); var++){
	title = vars[var].cuts; if(title=="1") title = "";
	title.ReplaceAll("nvmus==1&&nmus==1&&nvels==0","1 #mu");
	title.ReplaceAll("nvmus10==0&&nvels10==0", "0 leptons");  
	title.ReplaceAll("(nmus+nels)", "n_{lep}");  title.ReplaceAll("njets30","n_{jets}^{30}"); 
	title.ReplaceAll("els_pt","p^{e}_{T}");title.ReplaceAll("mus_pt","p^{#mu}_{T}");
	title.ReplaceAll("mus_reliso","RelIso"); title.ReplaceAll("els_reliso","RelIso");
	title.ReplaceAll("mus_miniso_tr15","MiniIso"); title.ReplaceAll("els_miniso_tr15","MiniIso");
	title.ReplaceAll("njets","n_{jets}");title.ReplaceAll("abs(lep_id)==13&&","");
	title.ReplaceAll(">=", " #geq "); title.ReplaceAll(">", " > "); title.ReplaceAll("&&", ", "); 
	title.ReplaceAll("met", "MET"); title.ReplaceAll("ht", "H_{T}");  title.ReplaceAll("mt", "m_{T}"); 
	title.ReplaceAll("nleps==1", "1 lepton");  title.ReplaceAll("nbm","n_{b}"); title.ReplaceAll("==", " = "); 
	title.ReplaceAll("nbl[1]","n_{b,l}");
	hname = "histo"; hname += var; hname += his; hname += sam;
	varhisto[his][sam].push_back(new TH1D(hname, title, vars[var].nbins, vars[var].minx, vars[var].maxx));
	varhisto[his][sam][var]->SetXTitle(vars[var].title);
	varhisto[his][sam][var]->SetYTitle(vars[var].varname);
	varhisto[his][sam][var]->SetLineColor(Samples[sam].color);
	varhisto[his][sam][var]->SetLineStyle(Samples[sam].style);
	varhisto[his][sam][var]->SetLineWidth(3);
	varhisto[his][sam][var]->SetMarkerStyle(20);
	varhisto[his][sam][var]->SetMarkerSize(1.7);
	varhisto[his][sam][var]->SetMarkerColor(Samples[sam].color);
	varhisto[his][sam][var]->Sumw2();
      } // Loop over variables to plot
    } // Loop over histo type: 1->MJ, 2->w, 3->w*MJ^2

    // Looping over entries
    small_tree tree(Samples[sam].file[0].Data());
    // // For now, will use only one sample per small_tree
    // for(unsigned insam(1); insam < Samples[vars[0].samples[sam]].file.size(); insam++)
    //   tree.Add(Samples[vars[0].samples[sam]].file[insam]);
    long nentries = tree.GetEntries();
    //nentries = 5;
    //cout<<endl<<"Doing "<<Samples[sam].file[0].Data()<<endl;
    for (long entry(0); entry<nentries; entry++){
      tree.GetEntry(entry);

      if(tree.ht()<500 || tree.met()<200 || tree.njets()<4 || (tree.nmus()+tree.nels())!=0) continue;
      for(unsigned var(0); var<vars.size(); var++){
	double weight(tree.weight()*4), ntrupv(tree.ntrupv()), mj(0);
	switch(var){
	case 0:
	  mj = tree.mj_10();
	  break;
	case 1:
	  mj = tree.mj_30();
	  break;
	case 2:
	  mj = tree.mj_40();
	  break;
	case 3:
	  mj = tree.mj_cands();
	  break;
	case 4:
	  mj = tree.mj_cands_trim();
	  break;
	default:
	  break;
	}
	varhisto[0][sam][var]->Fill(ntrupv, weight*mj);
	varhisto[1][sam][var]->Fill(ntrupv, weight);
	varhisto[2][sam][var]->Fill(ntrupv, weight*mj*mj);
	varhisto[3][sam][var]->Fill(ntrupv, weight*weight);
	//cout<<entry<<": mj "<<mj<<", w "<<weight<<", npv "<<ntrupv<<endl;
      } // Loop over variables to plot

    } // Loop over small_ntuple entries
  } // Loop over samples

  // Making the plots
  for(unsigned var(0); var<vars.size(); var++){
    for(unsigned sam(0); sam < Samples.size(); sam++){
      for(int bin(1); bin<=varhisto[0][sam][var]->GetNbinsX(); bin++){
	double mean(varhisto[0][sam][var]->GetBinContent(bin));
	double denom(varhisto[1][sam][var]->GetBinContent(bin));
	double wX2(varhisto[2][sam][var]->GetBinContent(bin));
	double w2(varhisto[3][sam][var]->GetBinContent(bin));
	double effn(pow(denom,2)); // Effective number of entries, see TH1::GetEffectiveEntries()
	if(denom) {
	  mean /= denom;
	  wX2 /= denom;
	  effn /= w2;
	} else {
	  mean = 0;
	  wX2 = 0;
	  effn = 1;
	}
	varhisto[0][sam][var]->SetBinContent(bin, mean);
	varhisto[0][sam][var]->SetBinError(bin, sqrt((wX2-mean*mean)/effn));	  
      } // Loop over number of bins

      varhisto[0][sam][var]->Fit(&linfit,"M Q N","", vars[var].minx, vars[var].maxx);
      leghisto = Samples[sam].label+" [<M_{J}> = "+RoundNumber(linfit.GetParameter(0),1)+" + "+
	RoundNumber(linfit.GetParameter(1),1)+"n_{PV}]";
      leg.AddEntry(varhisto[0][sam][var], leghisto);
      // cout<<Samples[sam].label<<" [<M_{J}> = ("<<RoundNumber(linfit.GetParameter(0),1)<<" +- "<<
      // 	RoundNumber(linfit.GetParError(0),1)<<") + ("<<
      // 	RoundNumber(linfit.GetParameter(1),1)<<" +- "<<RoundNumber(linfit.GetParError(1),1)<<")n_{PV}]"<<endl;
      if(sam==0){
	varhisto[0][sam][var]->SetMinimum(0);
	varhisto[0][sam][var]->SetMaximum(1650);
	varhisto[0][sam][var]->Draw("e");
      } else varhisto[0][sam][var]->Draw("e same");
    }  // Loop over samples
    leg.Draw();
    pname = "plots/pu_"+vars[var].tag+".eps";
    can.SaveAs(pname);  
    leg.Clear();
  } // Loop over variables to plot

  // Deleting pointers
  for(unsigned var(0); var<vars.size(); var++)
    for(unsigned sam(0); sam < Samples.size(); sam++)
      for(unsigned his(0); his < 4; his++)    
	if(varhisto[his][sam][var]) varhisto[his][sam][var]->Delete();
}

