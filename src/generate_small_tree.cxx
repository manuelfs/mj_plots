#include "generate_small_tree.hpp"

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

int main(){
  std::ofstream cppFile("src/small_tree.cpp"), hppFile("inc/small_tree.hpp");

  vector<string> variables;
  //////////////////   Provenance   /////////////
  variables.push_back("int run");
  variables.push_back("int lumiblock");
  variables.push_back("int event");

  ////////////////////   Global   ///////////////////
  variables.push_back("float weight");
  variables.push_back("int npv");
  variables.push_back("int ntrupv");
  variables.push_back("float ntrupv_mean");
  variables.push_back("float met");
  variables.push_back("float met_phi");
  variables.push_back("float ht");
  variables.push_back("float ht30");
  variables.push_back("float mht30");
  variables.push_back("float mt");
  variables.push_back("float dphi_wlep");
  variables.push_back("float mindphin_metjet");
  variables.push_back("float min_mt_bmet");
  variables.push_back("float min_mt_bmet_with_w_mass");

  //////////////////   Leptons   ////////////////
  variables.push_back("int nleps");
  variables.push_back("int nmus");
  variables.push_back("int nvmus");
  variables.push_back("int nvmus10");
  variables.push_back("std::vector<float> mus_pt");
  variables.push_back("std::vector<float> mus_eta");
  variables.push_back("std::vector<float> mus_phi");
  variables.push_back("std::vector<bool> mus_sigid");
  variables.push_back("std::vector<float> mus_reliso");
  variables.push_back("std::vector<float> mus_reliso_r04"); // recomputed default
  variables.push_back("std::vector<float> mus_reliso_r03"); // R_iso = 0.3
  variables.push_back("std::vector<float> mus_reliso_r02"); // R_iso = 0.2
  variables.push_back("std::vector<float> mus_reliso_r01"); // R_iso = 0.1
  variables.push_back("std::vector<float> mus_reliso_r015"); // R_iso = 0.15
  variables.push_back("std::vector<float> mus_miniso_10_ch"); // R_iso = 10./pT, charged PF only
  variables.push_back("std::vector<float> mus_miniso_15");   // R_iso = 15./pT
  variables.push_back("std::vector<float> mus_miniso_tr10"); // R_iso = max(0.05,min(0.3,10./pT))
  variables.push_back("std::vector<float> mus_miniso_tr07"); // R_iso = max(0.05,min(0.3,7./pT))
  variables.push_back("std::vector<float> mus_miniso_tr07_ch"); // R_iso = max(0.05,min(0.3,7./pT)), charged PF cands only
  variables.push_back("std::vector<float> mus_miniso_tr10_pfpu"); // R_iso = max(0.05,min(0.3,10./pT)), PF PU
  variables.push_back("std::vector<float> mus_ptrel_0");
  variables.push_back("std::vector<float> mus_mindr_0");
  variables.push_back("std::vector<float> mus_ptrel_25");
  variables.push_back("std::vector<float> mus_mindr_25");
  variables.push_back("std::vector<float> mus_ptrel_rem_0");
  variables.push_back("std::vector<float> mus_mindr_rem_0");
  variables.push_back("std::vector<float> mus_ptrel_rem_25");
  variables.push_back("std::vector<float> mus_mindr_rem_25");
  variables.push_back("std::vector<float> mus_tru_dr");
  variables.push_back("std::vector<int> mus_tru_id");
  variables.push_back("std::vector<int> mus_tru_momid");
  variables.push_back("std::vector<bool> mus_tru_tm");
  variables.push_back("std::vector<int> mus_charge");

  variables.push_back("int nels");
  variables.push_back("int nvels");
  variables.push_back("int nvels10");
  variables.push_back("std::vector<float> els_pt");
  variables.push_back("std::vector<float> els_eta");
  variables.push_back("std::vector<float> els_phi");
  variables.push_back("std::vector<bool> els_sigid");
  variables.push_back("std::vector<bool> els_ispf");
  variables.push_back("std::vector<float> els_reliso");
  variables.push_back("std::vector<float> els_reliso_r04"); // R_iso = 0.4
  variables.push_back("std::vector<float> els_reliso_r03"); // recomputed default
  variables.push_back("std::vector<float> els_reliso_r02"); // R_iso = 0.2
  variables.push_back("std::vector<float> els_reliso_r01"); // R_iso = 0.1
  variables.push_back("std::vector<float> els_reliso_r015"); // R_iso = 0.15
  variables.push_back("std::vector<float> els_miniso_10_ch"); // R_iso = 10./pT, charged PF only
  variables.push_back("std::vector<float> els_miniso_15");   // R_iso = 15./pT
  variables.push_back("std::vector<float> els_miniso_tr10"); // R_iso = max(0.05,min(0.3,10./pT))
  variables.push_back("std::vector<float> els_miniso_tr07"); // R_iso = max(0.05,min(0.3,7./pT))
  variables.push_back("std::vector<float> els_miniso_tr07_ch"); // R_iso = max(0.05,min(0.3,7./pT)), charged PF cands only
  variables.push_back("std::vector<float> els_miniso_tr10_pfpu"); // R_iso = max(0.05,min(0.3,10./pT)), PF PU
  variables.push_back("std::vector<float> els_ptrel_0");
  variables.push_back("std::vector<float> els_mindr_0");
  variables.push_back("std::vector<float> els_ptrel_25");
  variables.push_back("std::vector<float> els_mindr_25");
  variables.push_back("std::vector<float> els_ptrel_rem_0");
  variables.push_back("std::vector<float> els_mindr_rem_0");
  variables.push_back("std::vector<float> els_ptrel_rem_25");
  variables.push_back("std::vector<float> els_mindr_rem_25");
  variables.push_back("std::vector<float> els_tru_dr");
  variables.push_back("std::vector<int> els_tru_id");
  variables.push_back("std::vector<int> els_tru_momid");
  variables.push_back("std::vector<bool> els_tru_tm");
  variables.push_back("std::vector<int> els_charge");

  ///////////////////   Taus   //////////////////
  variables.push_back("int ntaus");
  variables.push_back("int ntaus_mt100");
  variables.push_back("int ntaus_againstEMu");
  variables.push_back("int ntaus_againstEMu_mt100");
  variables.push_back("std::vector<float> taus_pt");
  variables.push_back("std::vector<bool> taus_againstEMu");
  variables.push_back("std::vector<float> taus_eta");
  variables.push_back("std::vector<float> taus_phi");
  variables.push_back("std::vector<float> taus_chargedIsoPtSum");
  variables.push_back("std::vector<float> taus_mt");
  variables.push_back("std::vector<int> taus_n_pfcands");
  variables.push_back("std::vector<int> taus_decayMode");
  variables.push_back("std::vector<float> taus_CombinedIsolationDeltaBetaCorrRaw3Hits");

  ///////////////////   Tracks   //////////////////
  variables.push_back("std::vector<int> tks_id");
  variables.push_back("std::vector<float> tks_pt");
  variables.push_back("std::vector<float> tks_eta");
  variables.push_back("std::vector<float> tks_phi");
  variables.push_back("std::vector<float> tks_r03_ch");
  variables.push_back("std::vector<float> tks_r02_ch");
  variables.push_back("std::vector<float> tks_mini_ch");
  variables.push_back("std::vector<float> tks_r03_ne");
  variables.push_back("std::vector<float> tks_r02_ne");
  variables.push_back("std::vector<float> tks_mini_ne");
  variables.push_back("std::vector<int> tks_tru_id");
  variables.push_back("std::vector<float> tks_tru_dr");
  variables.push_back("std::vector<float> tks_tru_dp");
  variables.push_back("std::vector<bool> tks_from_w");
  variables.push_back("std::vector<bool> tks_from_tau");
  variables.push_back("std::vector<bool> tks_from_taulep");
  variables.push_back("std::vector<bool> tks_from_tauhad");
  variables.push_back("std::vector<int> tks_num_prongs");

  ///////////////////   Isolated tracks   //////////////////
  variables.push_back("int nisotk10");
  variables.push_back("int nisotk15");
  variables.push_back("int nisotk10_mt100");
  variables.push_back("int nisotk15_mt100");
  variables.push_back("std::vector<float> isotk_pt");
  variables.push_back("std::vector<float> isotk_iso");
  variables.push_back("std::vector<float> isotk_mt");
  variables.push_back("std::vector<float> isotk_eta");
  variables.push_back("std::vector<float> isotk_dzpv");

  ///////////////////   Jets   //////////////////
  variables.push_back("int njets");
  variables.push_back("int nbl");
  variables.push_back("int nbm");
  variables.push_back("int nbt");
  variables.push_back("int njets30");
  variables.push_back("int nbl30");
  variables.push_back("int nbm30");
  variables.push_back("int nbt30");

  variables.push_back("std::vector<float> jets_pt");
  variables.push_back("std::vector<float> jets_eta");
  variables.push_back("std::vector<float> jets_phi");
  variables.push_back("std::vector<float> jets_csv");
  variables.push_back("std::vector<float> jets_id");
  variables.push_back("std::vector<float> jets_dphi_lep");
  variables.push_back("std::vector<float> jets_dphi_met");
  variables.push_back("std::vector<float> jets_dphi_sum");
  variables.push_back("std::vector<bool> jets_bhad");
  variables.push_back("std::vector<bool> jets_blep");
  variables.push_back("std::vector<bool> jets_bhad_tru");
  variables.push_back("std::vector<bool> jets_blep_tru");

  ////////////////////   MC   ///////////////////
  variables.push_back("unsigned mc_type");
  variables.push_back("std::vector<float>  mc_pt");
  variables.push_back("std::vector<float>  mc_eta");
  variables.push_back("std::vector<float>  mc_phi");
  variables.push_back("std::vector<int>    mc_id");
  variables.push_back("std::vector<size_t> mc_mom");
  //variables.push_back("std::vector<float> mc_momid");
  //variables.push_back("std::vector<float> mc_gmomid");

  ///////////////////   Fat Jets   //////////////
  variables.push_back("int nfjets_10");  // for fjets_pt>50, all with skinny jets pt>10
  variables.push_back("float mj_10");    // for fjets_pt>50, all with skinny jets pt>10
  variables.push_back("std::vector<float> fjets_10_pt");
  variables.push_back("std::vector<float> fjets_10_eta");
  variables.push_back("std::vector<float> fjets_10_phi");
  variables.push_back("std::vector<float> fjets_10_m");

  variables.push_back("int nfjets_20");  // for fjets_pt>50, all with skinny jets pt>20
  variables.push_back("float mj_20");    // for fjets_pt>50, all with skinny jets pt>20
  variables.push_back("std::vector<float> fjets_20_pt");
  variables.push_back("std::vector<float> fjets_20_eta");
  variables.push_back("std::vector<float> fjets_20_phi");
  variables.push_back("std::vector<float> fjets_20_m");

  variables.push_back("int nfjets_30");  // for fjets_pt>50, all with skinny jets pt>30
  variables.push_back("float mj_30");    // for fjets_pt>50, all with skinny jets pt>30
  variables.push_back("std::vector<float> fjets_30_pt");
  variables.push_back("std::vector<float> fjets_30_eta");
  variables.push_back("std::vector<float> fjets_30_phi");
  variables.push_back("std::vector<float> fjets_30_m");

  variables.push_back("int nfjets_40");  // for fjets_pt>50, all with skinny jets pt>40
  variables.push_back("float mj_40");    // for fjets_pt>50, all with skinny jets pt>40
  variables.push_back("std::vector<float> fjets_40_pt");
  variables.push_back("std::vector<float> fjets_40_eta");
  variables.push_back("std::vector<float> fjets_40_phi");
  variables.push_back("std::vector<float> fjets_40_m");

  variables.push_back("int nfjets_nolep_30");  // for fjets_pt>50, all with skinny jets pt>30, no leptons
  variables.push_back("float mj_nolep_30");    // for fjets_pt>50, all with skinny jets pt>30, no leptons
  variables.push_back("std::vector<float> fjets_nolep_30_pt");
  variables.push_back("std::vector<float> fjets_nolep_30_eta");
  variables.push_back("std::vector<float> fjets_nolep_30_phi");
  variables.push_back("std::vector<float> fjets_nolep_30_m");

  variables.push_back("int nfjets_siglep_30");  // for fjets_pt>50, all with skinny jets pt>30, with sig leptons
  variables.push_back("float mj_siglep_30");    // for fjets_pt>50, all with skinny jets pt>30, with sig leptons
  variables.push_back("std::vector<float> fjets_siglep_30_pt");
  variables.push_back("std::vector<float> fjets_siglep_30_eta");
  variables.push_back("std::vector<float> fjets_siglep_30_phi");
  variables.push_back("std::vector<float> fjets_siglep_30_m");

  variables.push_back("int nfjets_cands");  // for fjets_pt>50, all pfcands
  variables.push_back("float mj_cands");    // for fjets_pt>50, all pfcands
  variables.push_back("std::vector<float> fjets_cands_pt");
  variables.push_back("std::vector<float> fjets_cands_eta");
  variables.push_back("std::vector<float> fjets_cands_phi");
  variables.push_back("std::vector<float> fjets_cands_m");

  variables.push_back("int nfjets_cands_trim");  // for fjets_pt>50, all pfcands
  variables.push_back("float mj_cands_trim");    // for fjets_pt>50, all pfcands
  variables.push_back("std::vector<float> fjets_cands_trim_pt");
  variables.push_back("std::vector<float> fjets_cands_trim_eta");
  variables.push_back("std::vector<float> fjets_cands_trim_phi");
  variables.push_back("std::vector<float> fjets_cands_trim_m");

  variables.push_back("int nfjets_r08");  // for fjets_pt>50, radius 0.8
  variables.push_back("float mj_r08");    // for fjets_pt>50, radius 0.8
  variables.push_back("std::vector<float> fjets_r08_pt");
  variables.push_back("std::vector<float> fjets_r08_eta");
  variables.push_back("std::vector<float> fjets_r08_phi");
  variables.push_back("std::vector<float> fjets_r08_m");

  variables.push_back("int nfjets_r10");  // for fjets_pt>50, radius 1.0
  variables.push_back("float mj_r10");    // for fjets_pt>50, radius 1.0
  variables.push_back("std::vector<float> fjets_r10_pt");
  variables.push_back("std::vector<float> fjets_r10_eta");
  variables.push_back("std::vector<float> fjets_r10_phi");
  variables.push_back("std::vector<float> fjets_r10_m");

  variables.push_back("int nfjets_r14");  // for fjets_pt>50, radius 1.4
  variables.push_back("float mj_r14");    // for fjets_pt>50, radius 1.4
  variables.push_back("std::vector<float> fjets_r14_pt");
  variables.push_back("std::vector<float> fjets_r14_eta");
  variables.push_back("std::vector<float> fjets_r14_phi");
  variables.push_back("std::vector<float> fjets_r14_m");

  variables.push_back("int nfjets_eta25");  // for fjets_pt>50, radius 1.5
  variables.push_back("float mj_eta25");    // for fjets_pt>50, radius 1.5
  variables.push_back("std::vector<float> fjets_eta25_pt");
  variables.push_back("std::vector<float> fjets_eta25_eta");
  variables.push_back("std::vector<float> fjets_eta25_phi");
  variables.push_back("std::vector<float> fjets_eta25_m");

  vector<Variable> vars = ParseVariables(variables);

  //////////  hpp file  ///////////////////
  hppFile << "// small_tree: class that contains reduced cfA ntuples\n";
  hppFile << "// File generated with generate_small_tree.exe\n\n";
  hppFile << "#ifndef H_SMALL_TREE\n";
  hppFile << "#define H_SMALL_TREE\n\n";

  hppFile << "#include <vector>\n";
  hppFile << "#include <string>\n\n";

  hppFile << "#include \"TTree.h\"\n";
  hppFile << "#include \"TChain.h\"\n\n";

  hppFile << "class small_tree{\n";
  hppFile << "public:\n";
  hppFile << "  small_tree(); // Constructor to create tree\n";
  hppFile << "  small_tree(const std::string &filename); // Constructor to read tree\n\n";

  hppFile << "  long GetEntries() const;\n";
  hppFile << "  void GetEntry(const long entry);\n\n";

  hppFile << "  void Fill();\n";
  hppFile << "  void Write();\n\n";

  hppFile << "  static const double bad_val_;\n\n";

  hppFile << "  ~small_tree();\n\n";

  for(vector<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    hppFile << "  " << var->type << " const & " << var->name << "() const;\n";
    hppFile << "  " << var->type << " & " << var->name << "();\n";
  }

  hppFile << "private:\n";
  hppFile << "  TChain chain_;\n";
  hppFile << "  TTree tree_;\n";
  hppFile << "  long entry_;\n";
  hppFile << "  const bool read_only_;\n\n";
  for(vector<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    hppFile << "  " << var->type << ' ' << var->name << "_;\n";
    if (var->type.find("vector<")!=string::npos) hppFile << "  " << var->type << " *p_" << var->name << "_;\n";
    hppFile << "  TBranch *b_" << var->name << "_;\n";
    hppFile << "  mutable bool c_" << var->name << "_;\n";
  }

  hppFile << "};\n\n";
  hppFile << "#endif" << endl;

  //////////  cpp file  ///////////////////
  cppFile << "// small_tree: class that contains reduced cfA ntuples\n";
  cppFile << "// File generated with generate_small_tree.exe\n\n";
  cppFile << "#include \"small_tree.hpp\"\n\n";

  cppFile << "#include <iostream>\n";
  cppFile << "#include <stdexcept>\n";
  cppFile << "#include <string>\n";
  cppFile << "#include <vector>\n\n";

  cppFile << "#include \"TTree.h\"\n";
  cppFile << "#include \"TBranch.h\"\n";
  cppFile << "#include \"TChain.h\"\n\n";

  cppFile << "using namespace std;\n\n";

  cppFile << "const double small_tree::bad_val_ = -999.;\n\n";

  cppFile << "small_tree::small_tree():\n";
  cppFile << "  chain_(\"junk\", \"junk\"),\n";
  cppFile << "  tree_(\"tree\", \"tree\"),\n";
  cppFile << "  entry_(0),\n";
  cppFile << "  read_only_(false),\n";
  for(size_t i = 0; i < vars.size()-1; ++i){
    cppFile << "  " << vars.at(i).name << "_(0),\n";
    cppFile << "  b_" << vars.at(i).name << "_(tree_.Branch(\"" <<vars.at(i).name << "\", &" << vars.at(i).name << "_)),\n";
    cppFile << "  c_" << vars.at(i).name << "_(false),\n";
  }
  cppFile << "  " << vars.back().name << "_(0),\n";
  cppFile << "  b_" << vars.back().name << "_(tree_.Branch(\"" <<vars.back().name << "\", &" << vars.back().name << "_)),\n";
  cppFile << "  c_" << vars.back().name << "_(false){\n";
  cppFile << "}\n\n";

  cppFile << "small_tree::small_tree(const string &filename):\n";
  cppFile << "  chain_(\"tree\",\"tree\"),\n";
  cppFile << "  tree_(\"junk\",\"junk\"),\n";
  cppFile << "  entry_(0),\n";
  cppFile << "  read_only_(true),\n";
  for(size_t i = 0; i < vars.size()-1; ++i){
    cppFile << "  " << vars.at(i).name << "_(0),\n";
    if (vars.at(i).type.find("vector<")!=string::npos) cppFile << "  p_" << vars.at(i).name << "_(&" << vars.at(i).name << "_),\n";
    cppFile << "  b_" << vars.at(i).name << "_(NULL),\n";
    cppFile << "  c_" << vars.at(i).name << "_(false),\n";
  }
  cppFile << "  " << vars.back().name << "_(0),\n";
  if (vars.back().type.find("vector<")!=string::npos) cppFile << "  p_" << vars.back().name << "_(&" << vars.back().name << "_),\n";
  cppFile << "  b_" << vars.back().name << "_(NULL),\n";
  cppFile << "  c_" << vars.back().name << "_(false){\n";
  cppFile << "  chain_.Add(filename.c_str());\n";
  for(vector<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    if (var->type.find("vector<")!=string::npos) cppFile << "  chain_.SetBranchAddress(\"" << var->name << "\", &p_" << var->name << "_, &b_" << var->name << "_);\n";
    else cppFile << "  chain_.SetBranchAddress(\"" << var->name << "\", &" << var->name << "_, &b_" << var->name << "_);\n";
  }
  cppFile << "}\n\n";

  cppFile << "void small_tree::Fill(){\n";
  cppFile << "  if(read_only_){\n";
  cppFile << "    throw std::logic_error(\"Trying to write to read-only tree\");\n";
  cppFile << "  }else{\n";
  cppFile << "    tree_.Fill();\n";
  cppFile << "  }\n";

  cppFile << "  //Resetting variables\n";
  for(vector<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    if(Contains(var->type, "vector")){
      cppFile << "  " << var->name << "_.clear();\n";
    }else if(Contains(var->type, "tring")){
      cppFile << "  " << var->name << "_ = \"\";\n";
    }else{
      cppFile << "  " << var->name << "_ = static_cast<" << var->type << ">(bad_val_);\n";
    }
  }
  cppFile << "}\n\n";

  cppFile << "void small_tree::Write(){\n";
  cppFile << "  if(read_only_){\n";
  cppFile << "    chain_.Write();\n";
  cppFile << "  }else{\n";
  cppFile << "    tree_.Write();\n";
  cppFile << "  }\n";
  cppFile << "}\n\n";

  cppFile << "small_tree::~small_tree(){\n";
  cppFile << "}\n\n";

  cppFile << "long small_tree::GetEntries() const{\n";
  cppFile << "  if(read_only_){\n";
  cppFile << "    return chain_.GetEntries();\n";
  cppFile << "  }else{\n";
  cppFile << "    return tree_.GetEntries();\n";
  cppFile << "  }\n";
  cppFile << "}\n\n";

  cppFile << "void small_tree::GetEntry(const long entry){\n";
  cppFile << "  if(!read_only_){\n";
  cppFile << "    throw std::logic_error(\"Trying to read from write-only tree.\");\n";
  cppFile << "  }\n\n";

  for(vector<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    cppFile << "  c_" << var->name << "_ = false;\n";
  }
  cppFile << "  entry_ = chain_.LoadTree(entry);\n";
  cppFile << "}\n";

  for(vector<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    cppFile << var->type << " const & small_tree::" << var->name << "() const{\n";
    cppFile << "  if(!read_only_){\n";
    cppFile << "    throw std::runtime_error(\"Trying to write to const tree.\");\n";
    cppFile << "  }\n";
    cppFile << "  if(!c_" << var->name << "_ && b_" << var->name <<"_){\n";
    cppFile << "    b_" << var->name << "_->GetEntry(entry_);\n";
    cppFile << "    c_" << var->name << "_ = true;\n";
    cppFile << "  }\n";
    cppFile << "  return " << var->name << "_;\n";
    cppFile << "}\n\n";
  }

  for(vector<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    cppFile << var->type << " & small_tree::" << var->name << "(){\n";
    cppFile << "  if(read_only_ && !c_" << var->name << "_ && b_" << var->name <<"_){\n";
    cppFile << "    b_" << var->name << "_->GetEntry(entry_);\n";
    cppFile << "    c_" << var->name << "_ = true;\n";
    cppFile << "  }\n";
    cppFile << "  return " << var->name << "_;\n";
    cppFile << "}\n\n";
  }

  cout<<"Written src/small_tree.cpp and inc/small_tree.hpp"<<endl;

  cppFile.close();
  hppFile.close();
  return 0;
}

vector<Variable> ParseVariables(const vector<string> &variables){
  vector<Variable> out(variables.size());
  for(size_t i = 0; i < variables.size(); ++i){
    string::size_type pos = variables.at(i).rfind(' ');
    if(pos != string::npos){
      out.at(i).type = variables.at(i).substr(0,pos);
      out.at(i).name = variables.at(i).substr(pos+1);
    }else{
      out.at(i).type = "";
      out.at(i).name = "";
      throw runtime_error("Bad variable: "+variables.at(i));
    }
  }
  return out;
}

bool Contains(const string &text, const string &pattern){
  return text.find(pattern) != string::npos;
}
