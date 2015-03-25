// plot_mj_correl: Macro that plots the correlation matrices coming out of TMVA

#define INT_ROOT
#include "inc/styles.hpp"
#include "src/styles.cpp"
#include "inc/utilities.hpp"
#include "src/utilities.cpp"

#include <vector> 
#include "TChain.h"
#include "TFile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"
#include "TColor.h"

using namespace std;
using std::cout;
using std::endl;

class sample_class {
public:
  sample_class(TString file, TString name, TString label, TString type, int color);
  TString file, name, label, type;
  int color;
};

sample_class::sample_class(TString ifile, TString iname, TString ilabel, TString itype, int icolor):
  file(ifile),
  name(iname),
  label(ilabel),
  type(itype),
  color(icolor){
}

void plot_mj_correl() { 
  styles style("2Dtitle"); style.setDefaultStyle();
  TCanvas can;
  TString folder="root/", pname;

  vector<sample_class> files;
  files.push_back(sample_class(folder+"bdt_T1tttt1500_TTJet_mj_30htnjets30met.root",
			       "T1tttt(1500,100)", "t1t","S",2));
  files.push_back(sample_class(folder+"bdt_T1tttt1200_TTJet_mj_30htnjets30met.root",
			       "T1tttt(1200,800)", "t1tc","S",3));
  files.push_back(sample_class(folder+"bdt_T1tttt1500_TTJet_mj_30htnjets30met.root",
			       "t#bar{t}", "tt","B",1));
  files.push_back(sample_class(folder+"bdt_T1tttt1500_QCD_Pt_mj_30htnjets30met.root",
			       "QCD", "qcd","B",4));

  const unsigned num_conts = 50;
  const unsigned num_stops = 2;
  double stops[num_stops] = {0., 1.};
  double red[num_stops] = {1., 0.3};
  double green[num_stops] = {1., 0.3};
  double blue[num_stops] = {1., 1.};

  for(unsigned ifile(0); ifile < files.size(); ifile++){
    TFile file(files[ifile].file);
    if (file.IsZombie()) {
      cout << "Error opening "<<files[ifile].file << "  -  Exiting "<<endl;
      continue;
    }
    file.cd();

    switch(files[ifile].color){
    case 1:
      red[1] = 0.39;
      green[1] = 0.51;
      blue[1] = 0.78;
      break;
    case 2:
      red[1] = 0.23;
      green[1] = 0.72;
      blue[1] = 0.52;
      break;
    case 3:
      red[1] = 0.78;
      green[1] = 0.47;
      blue[1] = 0.63;
      break;
    case 4:
      red[1] = 0.08;
      green[1] = 0.87;
      blue[1] = 0.83;
      break;
    default:
      break;
    }
    gStyle->SetNumberContours(num_conts);
    TColor::CreateGradientColorTable(num_stops, stops, red, green, blue, num_conts);

    TH2F* hcorrel = (TH2F*)file.Get("CorrelationMatrix"+files[ifile].type);
    hcorrel->SetMinimum(0);
    hcorrel->SetMaximum(100);
    hcorrel->SetTitle("Correlation matrix for "+files[ifile].name);
    hcorrel->SetZTitle("#rho (%)");
    hcorrel->SetLabelFont(style.nFont,"xyz");
    hcorrel->SetTitleFont(style.nFont,"xyz");
    hcorrel->SetTitleFont(style.nFont);
    hcorrel->SetTitleSize(style.TextSize);
    hcorrel->SetTitleSize(style.TextSize,"xyz");
    hcorrel->SetLabelSize(style.LabelSize*2,"xy");
    hcorrel->SetLabelSize(style.LabelSize*0.81,"z");
    hcorrel->SetMarkerSize(3);
    hcorrel->SetMarkerColor(1);
    hcorrel->SetTitleOffset(0.53,"z");

    hcorrel->Draw("colz");
    hcorrel->Draw("text same");
    
    pname = "plots/correlation_"+files[ifile].label+".pdf";
    can.SaveAs(pname);
    file.Close();
  }
  
}

