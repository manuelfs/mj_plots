// train_mj_bdt: trains various BDTs which combine MJ, njets, and HT. 

#include <stdexcept>
#include <iostream>
#include <vector>
#include <utility> 

#include "TFile.h"
#include "TString.h"
#include "TChain.h"
#include "TCut.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"

#include "utilities.hpp"

using namespace std;

class bdt_class {
public:
  bdt_class(vector<pair<TString, char> > ivars, TString isignal, TString ibkg);
  vector<pair<TString, char> > vars;
  TString signal, bkg, name;
};

bdt_class::bdt_class(vector<pair<TString, char> > ivars, TString isignal, TString ibkg):
  vars(ivars),
  signal(isignal),
  bkg(ibkg){
  name = "bdt_"+signal+"_"+bkg+"_";
  for(unsigned ind(0); ind < vars.size(); ind++) name += vars[ind].first;
  name.ReplaceAll("*","");
}

int main(){

  TString nTrain("20000"), nTest("0");
  TString ntuplefolder("archive/15-01-30/skim/"), rootfolder("root/");
  gSystem->mkdir(rootfolder, kTRUE);

  vector<pair<TString, char> > v_htmj;
  v_htmj.push_back(make_pair("ht",'F'));
  v_htmj.push_back(make_pair("mj_30",'F'));

  vector<pair<TString, char> > v_htnjets;
  v_htnjets.push_back(make_pair("ht",'F'));
  v_htnjets.push_back(make_pair("njets30",'I'));

  vector<pair<TString, char> > v_njetsmj;
  v_njetsmj.push_back(make_pair("njets30",'I'));
  v_njetsmj.push_back(make_pair("mj_30",'F'));

  vector<TString> v_signal;
  //v_signal.push_back("*T1tttt*1500*");
   v_signal.push_back("*T1tttt*1200*");
   v_signal.push_back("*T1qqqq*1400*");
   v_signal.push_back("*T1qqqq*1000*");
		      
  vector<TString> v_bkg;
  //v_bkg.push_back("*QCD_Pt*");
  v_bkg.push_back("*TTJet*");
  //v_bkg.push_back("*ZJetsToNuNu*");

  vector<bdt_class> bdts;
  for(unsigned isig(0); isig < v_signal.size(); isig++){
    for(unsigned ibkg(0); ibkg < v_bkg.size(); ibkg++){
      bdts.push_back(bdt_class(v_htmj, v_signal[isig], v_bkg[ibkg]));
      bdts.push_back(bdt_class(v_htnjets, v_signal[isig], v_bkg[ibkg]));
      bdts.push_back(bdt_class(v_njetsmj, v_signal[isig], v_bkg[ibkg]));
    }
  }

  for(unsigned ibdt(0); ibdt < bdts.size(); ibdt++){
    // Loading root file and trees
    TFile tmvaFile(rootfolder+bdts[ibdt].name+".root", "RECREATE");
    tmvaFile.cd();
    TChain signal("tree"), ttbar("tree");
    signal.Add(ntuplefolder+bdts[ibdt].signal);
    ttbar.Add(ntuplefolder+bdts[ibdt].bkg);

    // Creating TMVA factory
    TMVA::Factory *factory = new TMVA::Factory(bdts[ibdt].name, &tmvaFile,"!V:Silent:Color");
    for(unsigned ind(0); ind < bdts[ibdt].vars.size(); ind++) 
      factory->AddVariable(bdts[ibdt].vars[ind].first, bdts[ibdt].vars[ind].second);
    factory->AddSpectator("met := met", "MET", "GeV", 'F');
    factory->AddSpectator("ht := ht", "H_{T}", "GeV", 'F');
    factory->AddSpectator("mj_30 := mj_30", "M_{J}", "GeV", 'F');
    factory->AddSpectator("nvmus10 := nvmus10", "Number of veto muons", "", 'I');
    factory->AddSpectator("nvels10 := nvels10", "Number of veto electrons", "", 'I');
    factory->AddSpectator("njets30 := njets30", "Number of 30 GeV jets", "", 'I');
    factory->AddSignalTree    (&signal, 1.);
    factory->AddBackgroundTree(&ttbar, 1.);
    factory->SetSignalWeightExpression    ("weight");
    factory->SetBackgroundWeightExpression("weight");

    TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
    TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";
    factory->PrepareTrainingAndTestTree( mycuts, mycutb, "nTrain_Signal="+nTrain+":nTrain_Background="+nTrain+
					 ":nTest_Signal="+nTest+":nTest_Background="+nTest+
					 ":SplitMode=Random:NormMode=NumEvents:!V" );
    factory->BookMethod( TMVA::Types::kBDT, "BDT",
			 "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );
    
    factory->TrainAllMethods();    // Train MVAs using the set of training events
    factory->TestAllMethods();     // Evaluate all MVAs using the set of test events
    factory->EvaluateAllMethods(); // Evaluate and compare performance of all configured MVAs

    delete factory;
    cout<<"\t\tWritten "<<rootfolder+bdts[ibdt].name+".root"<<endl;

  } // Loop over BDTs

}

