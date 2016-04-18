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

class var_class {
public:
  var_class(TString variable, TString name, char type, TString unit="GeV");
  TString variable, name;
  char type;
  TString unit;
};

var_class::var_class(TString ivariable, TString iname, char itype, TString iunit):
  variable(ivariable),
  name(iname),
  type(itype),
  unit(iunit){
}

class bdt_class {
public:
  bdt_class(vector<var_class> ivars, TString isignal, TString ibkg);
  vector<var_class> vars;
  TString signal, bkg, name;
};

bdt_class::bdt_class(vector<var_class> ivars, TString isignal, TString ibkg):
  vars(ivars),
  signal(isignal),
  bkg(ibkg){
  name = "bdt_"+signal+"_"+bkg+"_";
  for(unsigned ind(0); ind < vars.size(); ind++) name += vars[ind].variable;
  name.ReplaceAll("*","");
}

int main(){

  TString nTrain("20000"), nTest("0");
  TString bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))  
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder
  
  TString ntuplefolder=bfolder+"/cms2r0/babymaker/babies/2015_01_30/small_tree/skim_ht500met200/";

  TString rootfolder(bfolder+"/cms2r0/babymaker/babies/2015_01_30/small_tree/bdt/");
  gSystem->mkdir(rootfolder, kTRUE);

  vector<var_class> v_htmj;
  v_htmj.push_back(var_class("ht","H_{T}",'F'));
  v_htmj.push_back(var_class("mj_30","M_{J}",'F'));

  vector<var_class> v_htnjets;
  v_htnjets.push_back(var_class("ht","H_{T}",'F'));
  v_htnjets.push_back(var_class("njets30","n_{jets}^{30}",'I'));

  vector<var_class> v_njetsmj;
  v_njetsmj.push_back(var_class("njets30","n_{jets}^{30}",'I'));
  v_njetsmj.push_back(var_class("mj_30","M_{J}",'F'));

  vector<var_class> v_njetsmjhtmet;
  v_njetsmjhtmet.push_back(var_class("mj_30","M_{J}",'F'));
  v_njetsmjhtmet.push_back(var_class("ht","H_{T}",'F'));
  v_njetsmjhtmet.push_back(var_class("njets30","n_{jets}^{30}",'I'));
  v_njetsmjhtmet.push_back(var_class("met","MET",'F'));

  vector<TString> v_signal;
  v_signal.push_back("*T1tttt*1500*");
  v_signal.push_back("*T1tttt*1200*");
  // v_signal.push_back("*T1qqqq*1400*");
  // v_signal.push_back("*T1qqqq*1000*");
		      
  vector<TString> v_bkg;
  // v_bkg.push_back("*QCD_Pt*");
  v_bkg.push_back("*TTJet*");
  //v_bkg.push_back("*ZJetsToNuNu*");

  vector<bdt_class> bdts;
  for(unsigned isig(0); isig < v_signal.size(); isig++){
    for(unsigned ibkg(0); ibkg < v_bkg.size(); ibkg++){
      bdts.push_back(bdt_class(v_htmj, v_signal[isig], v_bkg[ibkg]));
      bdts.push_back(bdt_class(v_htnjets, v_signal[isig], v_bkg[ibkg]));
      bdts.push_back(bdt_class(v_njetsmj, v_signal[isig], v_bkg[ibkg]));
      //bdts.push_back(bdt_class(v_njetsmjhtmet, v_signal[isig], v_bkg[ibkg]));
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
      factory->AddVariable(bdts[ibdt].vars[ind].variable, bdts[ibdt].vars[ind].name, 
			   bdts[ibdt].vars[ind].unit, bdts[ibdt].vars[ind].type);
    factory->AddSpectator("met := met", "MET", "GeV", 'F');
    factory->AddSpectator("ht := ht", "H_{T}", "GeV", 'F');
    factory->AddSpectator("mj_30 := mj_30", "M_{J}", "GeV", 'F');
    factory->AddSpectator("nvmus10 := nvmus10", "Number of veto muons", "", 'I');
    factory->AddSpectator("nvels10 := nvels10", "Number of veto electrons", "", 'I');
    factory->AddSpectator("nmus := nmus", "Number of muons", "", 'I');
    factory->AddSpectator("nels := nels", "Number of electrons", "", 'I');
    factory->AddSpectator("njets30 := njets30", "Number of 30 GeV jets", "", 'I');
    factory->AddSignalTree    (&signal, 1.);
    factory->AddBackgroundTree(&ttbar, 1.);
    factory->SetSignalWeightExpression    ("weight");
    factory->SetBackgroundWeightExpression("weight");

    TCut mycuts(""), mycutb("");
    factory->PrepareTrainingAndTestTree( mycuts, mycutb, "nTrain_Signal="+nTrain+":nTrain_Background="+nTrain+
					 ":nTest_Signal="+nTest+":nTest_Background="+nTest+
					 ":SplitMode=Random:NormMode=NumEvents:!V" );
    factory->BookMethod( TMVA::Types::kBDT, "BDT",
			 "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );
    
    factory->TrainAllMethods();    // Train MVAs using the set of training events
    factory->TestAllMethods();     // Evaluate all MVAs using the set of test events
    factory->EvaluateAllMethods(); // Evaluate and compare performance of all configured MVAs

    delete factory;
    tmvaFile.Write();
    tmvaFile.Close();
    cout<<"\t\t\t\t\t\t\t\t\t\t\tWritten "<<rootfolder+bdts[ibdt].name+".root"<<endl;

  } // Loop over BDTs

}

