// plot_roc: Macro that plots ROC curves
#ifndef H_PLOT_MJ_ROC
#define H_PLOT_MJ_ROC

#include <vector>
#include "TH1D.h"
#include "TString.h"

class marker_class {
public:
  marker_class(float icut, float isize, int icolor, int istyle);
  float cut, size;
  int color, style;
};

class var_class {
public:
  var_class(TString ivarname, float iminx, float imaxx, TString ititle, int icolor, 
	    int istyle=1, std::vector<marker_class> icuts=std::vector<marker_class>());
  TString title, varname;
  float minx, maxx;
  std::vector<marker_class> cuts;
  int color, style;
};

class sample_class {
public:
  sample_class(TString ilabel, std::vector<TString> ifiles, TString icut="1");
  TString label, cut;
  std::vector<TString> files;
};

TGraph MakeROC(TH1D &good, TH1D &bad, const bool less_is_better, std::vector<marker_class> cuts);
void DrawROC(std::vector<sample_class> samples, std::vector<var_class> vars, TString cuts, TString tag);



#endif
