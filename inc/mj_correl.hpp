#ifndef H_MJ_CORREL
#define H_MJ_CORREL

#include "small_tree_quick.hpp"
#include <string>
#include <vector>

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TString.h"
#include "TVector3.h"
#include "TLorentzVector.h"


class var{
  public:
    var(){};
    ~var(){};

    var(TString name_, TString axisttl_, int nbins_, double min_, double max_):
    name(name_),
    axisttl(axisttl_),
    nbins(nbins_),
    min(min_),
    max(max_){};

    TString name; 
    TString axisttl;
    int nbins;
    double min;
    double max;
};

class seln{
  public:
    seln(){};
    ~seln(){};

    seln(TString name_, 
         int nleps_, int nisotks_max_,
         double sig_lep_pt_min_, double veto_lep_pt_min_,
         double ht_min_,double ht_max_,
         double mj_min_,double mj_max_,
         double met_min_, double met_max_,
         double mt_min_, double mt_max_,
         int njets_min_, int njets_max_,
         int nbl_min_, int nbl_max_):
    name(name_), 
    nleps(nleps_), nisotks_max(nisotks_max_),
    sig_lep_pt_min(sig_lep_pt_min_),veto_lep_pt_min(veto_lep_pt_min_),
    ht_min(ht_min_),ht_max(ht_max_),
    mj_min(mj_min_),mj_max(mj_max_),
    met_min(met_min_), met_max(met_max_),
    mt_min(mt_min_), mt_max(mt_max_),
    njets_min(njets_min_), njets_max(njets_max_),
    nbl_min(nbl_min_), nbl_max(nbl_max_){};

    TString name;

    size_t nleps;
    size_t nisotks_max;

    double sig_lep_pt_min;
    double veto_lep_pt_min;
    double ht_min;
    double ht_max;
    double mj_min;
    double mj_max;
    double met_min;
    double met_max;
    double mt_min;
    double mt_max;

    int njets_min;
    int njets_max;
    int nbl_min;
    int nbl_max;
};

class sample{
  public:
    sample(){};
    ~sample(){};

    sample(TString filestr_, TString name_, bool has_tops_):
    filestr(filestr_),
    name(name_),
    has_tops(has_tops_){};

    TString filestr;
    TString name;
    bool has_tops;
};

class h1d{
  public:
    h1d(){};
    ~h1d(){};

    h1d(const var &var_, TString hname){
      hist = TH1D(hname, hname, var_.nbins, var_.min, var_.max); 
      hist.GetXaxis()->SetTitle(var_.axisttl);
      hist.GetYaxis()->SetTitle("Events");
      hist.SetLineWidth(2);
      hist.Sumw2();
    }

    TH1D hist;
    
    void fill(double value, double weight){
      double max_val = hist.GetBinCenter(hist.GetNbinsX());
      if (value<max_val) hist.Fill(value, weight);
      else hist.Fill(max_val, weight);
    };
};

class h2d{
  public:
    h2d(){};
    ~h2d(){};

    h2d(const var &xvar_, const var &yvar_, TString hname){
      hist = TH2D(hname, hname, 
               xvar_.nbins, xvar_.min, xvar_.max,
               yvar_.nbins, yvar_.min, yvar_.max); 
      hist.GetXaxis()->SetTitle(xvar_.axisttl);
      hist.GetYaxis()->SetTitle(yvar_.axisttl);
      // hist.SetMarkerStyle(20);
      // hist.SetMarkerSize(0.5);
      hist.Sumw2();
    }

    TH2D hist;

    void fill(double xval, double yval, double weight){
      hist.Fill(xval, yval, weight);
    };
};

class h3d{
  public:
    h3d(){};
    ~h3d(){};

    h3d(const var &xvar_, const var &yvar_, const var &zvar_, TString hname){
      hist = TH3D(hname, hname, 
               xvar_.nbins, xvar_.min, xvar_.max,
               yvar_.nbins, yvar_.min, yvar_.max, 
               zvar_.nbins, zvar_.min, zvar_.max); 
      hist.GetXaxis()->SetTitle(xvar_.axisttl);
      hist.GetYaxis()->SetTitle(yvar_.axisttl);
      hist.GetZaxis()->SetTitle(zvar_.axisttl);
      // hist.SetMarkerStyle(20);
      // hist.SetMarkerSize(0.5);
      hist.Sumw2();
    }

    TH3D hist;

    void fill(double xval, double yval, double zval, double weight){
      hist.Fill(xval, yval, zval, weight);
    };
};

void msgsvc(const unsigned &lvl, const TString &mymsg);

bool compMass(const TLorentzVector &fj1, const TLorentzVector &fj2);
bool compPT(const TVector3 &v1, const TVector3 &v2);
double getHT(std::vector<TVector3> &v);

bool passIsolation(const small_tree &tree, int ilep, bool isElectron, bool isveto, TString isotype, const double coneiso_cut);
bool passSelection(const small_tree &tree, const seln &iseln);

#endif
