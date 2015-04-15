#include "generate_small_tree.hpp"

#include <cstring>

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <set>

using namespace std;

string ToCaps(string str){
  for(string::iterator it = str.begin();
      it != str.end();
      ++it){
    *it = toupper(*it);
  }
  return str;
}

string execute(const string &cmd){
  FILE *pipe = popen(cmd.c_str(), "r");
  if(!pipe) throw runtime_error("Could not open pipe.");
  const size_t buffer_size = 128;
  char buffer[buffer_size];
  string result = "";
  while(!feof(pipe)){
    if(fgets(buffer, buffer_size, pipe) != NULL) result += buffer;
  }

  pclose(pipe);
  return result;
}

vector<string> Tokenize(const string& input,
                        const string& tokens = " "){
  char* ipt(new char[input.size()+1]);
  memcpy(ipt, input.data(), input.size());
  ipt[input.size()]=static_cast<char>(0);
  char* ptr(strtok(ipt, tokens.c_str()));
  vector<string> output(0);
  while(ptr!=NULL){
    output.push_back(ptr);
    ptr=strtok(NULL, tokens.c_str());
  }
  return output;
}

string FixName(string name){
  //Variable names can have alphanumeric characters and underscores only
  string allowed = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890_";

  //Remove illegal characters
  size_t pos = name.find_first_not_of(allowed);
  while(pos < name.size()){
    name.at(pos) = '_';
    pos = name.find_first_not_of(allowed);
  }

  //Replace double underscore with single underscore
  pos = name.rfind("__");
  while(pos < name.size()){
    name.replace(pos, 2, "_");
    pos = name.rfind("__");
  }

  //Remove leading and trailing spaces
  pos = 0;
  for(pos = 0; pos < name.size(); ++pos){
    if(name.at(pos) != ' ') break;
  }
  size_t endpos = name.size();
  for(endpos = name.size(); endpos != 0; --endpos){
    if(name.at(endpos-1) != ' ') break;
  }

  return name.substr(pos, endpos-pos);
}

set<Variable> GetVariables(const string &file_name){
  string allowed = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890_";
  set<Variable> vars;

  ifstream infile(("txt/small_tree_cfg/"+file_name).c_str());
  string line;
  while(getline(infile, line)){
    size_t start = line.find_first_not_of(" ");
    if(start >= line.size() || line.at(start) == '#' || line.at(start) == '/') continue;

    //Replace double space with single space
    size_t pos = line.rfind("  ");
    while(pos < line.size()){
      line.replace(pos, 2, " ");
      pos = line.rfind("  ");
    }
    size_t end = line.find_last_of(allowed)+1;
    size_t split = line.rfind(' ', end)+1;

    vars.insert(Variable(line.substr(start, split-start),
                         line.substr(split, end-split)));
  }
  infile.close();

  return vars;
}

int main(int argc, char *argv[]){
  int c = 0;
  bool do_event_handler = true;
  while((c=getopt(argc, argv, "t"))!=-1){
    switch(c){
    case 't':
      do_event_handler = false;
      break;
    default:
      break;
    }
  }
  vector<string> file_names = Tokenize(execute("ls txt/small_tree_cfg/ 2> /dev/null"), "\n");

  vector<pair<string, set<Variable> > > sep_vars(file_names.size());
  vector<string> fixed_names(file_names.size());
  for(size_t ifile = 0; ifile < file_names.size(); ++ifile){
    fixed_names.at(ifile) = FixName(file_names.at(ifile));
    sep_vars.at(ifile).first = fixed_names.at(ifile);
    sep_vars.at(ifile).second = GetVariables(file_names.at(ifile));
  }

  set<Variable> all_vars;
  for(size_t ifile = 0; ifile < sep_vars.size(); ++ifile){
    for(set<Variable>::const_iterator var = sep_vars.at(ifile).second.begin();
        var != sep_vars.at(ifile).second.end();
        ++var){
      all_vars.insert(*var);
    }
  }

  set<Variable> com_vars;
  if(sep_vars.size()){
    for(set<Variable>::const_iterator var = sep_vars.at(0).second.begin();
        var != sep_vars.at(0).second.end();
        ++var){
      bool found_in_all = true;
      for(size_t ifile = 1; found_in_all && ifile < sep_vars.size(); ++ifile){
        if(sep_vars.at(ifile).second.find(*var) == sep_vars.at(ifile).second.end()){
          found_in_all = false;
        }
      }
      if(found_in_all){
        com_vars.insert(*var);
      }
    }
    for(set<Variable>::const_iterator var = com_vars.begin();
        var != com_vars.end();
        ++var){
      for(size_t ifile = 0; ifile < sep_vars.size(); ++ifile){
        sep_vars.at(ifile).second.erase(*var);
      }
    }
  }

  WriteBaseHeader(all_vars, com_vars, fixed_names);
  WriteBaseSource(all_vars, com_vars, fixed_names);

  for(size_t ifile = 0; ifile < sep_vars.size(); ++ifile){
    WriteSepHeader(sep_vars.at(ifile));
    WriteSepSource(sep_vars.at(ifile));
  }

  if(do_event_handler){
    vector<string> existing = Tokenize(execute("ls src/event_handler* 2> /dev/null"));
    set<string> to_generate;
    for(size_t iname = 0; iname < fixed_names.size(); ++iname){
      bool found = false;
      for(size_t iexist = 0; !found && iexist < existing.size(); ++iexist){
        if(Contains(existing.at(iexist), fixed_names.at(iname))) found = true;
      }
      if(!found) to_generate.insert(fixed_names.at(iname));
    }

    GenerateEventHandlerBaseHeader();
    GenerateEventHandlerBaseSource();

    for(set<string>::const_iterator name = to_generate.begin();
        name != to_generate.end();
        ++name){
      GenerateEventHandlerHeader(*name);
      GenerateEventHandlerSource(*name);
    }

    GenerateEventHandlerMergeHeader(fixed_names);
    GenerateEventHandlerMergeSource(fixed_names);
  }
}

bool Contains(const string &text, const string &pattern){
  return text.find(pattern) != string::npos;
}

void WriteBaseHeader(const set<Variable> &all_vars,
                     const set<Variable> &com_vars,
                     const vector<string> &names){
  ofstream file("inc/small_tree.hpp");

  file << "// small_tree: base class to handle reduced tree ntuples\n";
  file << "// File generated with generate_small_tree.exe\n\n";

  file << "#ifndef H_SMALL_TREE\n";
  file << "#define H_SMALL_TREE\n\n";

  file << "#include <vector>\n";
  file << "#include <string>\n\n";

  file << "#include \"TTree.h\"\n";
  file << "#include \"TChain.h\"\n\n";

  file << "class small_tree{\n";
  file << "public:\n";
  file << "  small_tree(); // Constructor to create tree\n";
  file << "  small_tree(const std::string &filename); // Constructor to read tree\n\n";

  file << "  int Add(const std::string &filename);\n";

  file << "  long GetEntries() const;\n";
  file << "  virtual void GetEntry(const long entry);\n\n";

  file << "  virtual void Fill();\n";
  file << "  void Write();\n\n";

  file << "  virtual std::string Type() const;\n\n";

  file << "  static const double bad_val_;\n\n";

  file << "  virtual ~small_tree();\n\n";

  for(set<Variable>::const_iterator var = com_vars.begin();
      var != com_vars.end();
      ++var){
    file << "  " << var->type_ << " const & " << var->name_ << "() const;\n";
    file << "  " << var->type_ << " & " << var->name_ << "();\n";
  }
  file << '\n';

  for(set<Variable>::const_iterator var = all_vars.begin();
      var != all_vars.end();
      ++var){
    if(com_vars.find(*var) != com_vars.end()) continue;
    file << "  __attribute__((noreturn)) virtual "
         << var->type_ << " const & " << var->name_ << "() const;\n";
    file << "  __attribute__((noreturn)) virtual "
         << var->type_ << " & " << var->name_ << "();\n";
  }
  file << '\n';

  file << "protected:\n";
  file << "  TChain chain_;\n";
  file << "  TTree tree_;\n";
  file << "  long entry_;\n";
  file << "  const bool read_only_;\n\n";

  file << "private:\n";
  file << "  class VectorLoader{\n";
  file << "  public:\n";
  file << "    VectorLoader();\n";
  file << "  private:\n";
  file << "    static bool loaded_;\n";
  file << "  };\n\n";

  file << "  static VectorLoader vl_;\n";
  for(set<Variable>::const_iterator var = com_vars.begin();
      var != com_vars.end();
      ++var){
    file << "  " << var->type_ << ' ' << var->name_ << "_;\n";
    if(Contains(var->type_, "vector")){
      file << "  " << var->type_ << " *p_" << var->name_ << "_;\n";
    }
    file << "  TBranch *b_" << var->name_ << "_;\n";
    file << "  mutable bool c_" << var->name_ << "_;\n";
  }
  file << "};\n\n";

  file << "small_tree* NewTree(const std::type_info &type);\n\n";

  for(size_t i = 0; i < names.size(); ++i){
    file << "#include \"small_tree_" << names.at(i) << ".hpp\"\n";
  }
  file <<'\n';

  file << "#endif" << endl;

  file.close();
}

void WriteBaseSource(const set<Variable> &all_vars,
                     const set<Variable> &com_vars,
                     const vector<string> &names){
  ofstream file("src/small_tree.cpp");

  file << "// small_tree: base class to handle reduce tree ntuples\n";
  file << "//File generated with generate_small_tree.exe\n\n";

  file << "#include \"small_tree.hpp\"\n\n";

  file << "#include <stdexcept>\n";
  file << "#include <string>\n";
  file << "#include <iostream>\n";
  file << "#include <vector>\n\n";

  file << "#include \"TROOT.h\"\n";
  file << "#include \"TTree.h\"\n";
  file << "#include \"TBranch.h\"\n";
  file << "#include \"TChain.h\"\n\n";

  file << "using namespace std;\n\n";

  file << "bool small_tree::VectorLoader::loaded_ = false;\n\n";

  file << "small_tree::VectorLoader small_tree::vl_ = small_tree::VectorLoader();\n\n";

  file << "small_tree::VectorLoader::VectorLoader(){\n";
  file << "  if(!loaded_){\n";
  file << "    gROOT->ProcessLine(\"#include <vector>\");\n";
  file << "    loaded_ = true;\n";
  file << "  }\n";
  file << "}\n\n";

  file << "const double small_tree::bad_val_ = -999.;\n\n";

  file << "small_tree::small_tree():\n";
  file << "  chain_(\"junk\", \"junk\"),\n";
  file << "  tree_(\"tree\", \"tree\"),\n";
  file << "  entry_(0),\n";
  if(com_vars.size()){
    const set<Variable>::const_iterator com_end_2 = --com_vars.end();
    file << "  read_only_(false),\n";
    for(set<Variable>::const_iterator var = com_vars.begin();
        var != com_end_2;
        ++var){
      file << "  " << var->name_ << "_(0),\n";
      if(Contains(var->type_, "vector")){
        file << "  p_" << var->name_ << "_(&" << var->name_ << "_),\n";
        file << "  b_" << var->name_ << "_(tree_.Branch(\"" << var->name_ << "\", &p_" << var->name_ << "_)),\n";
      }else{
        file << "  b_" << var->name_ << "_(tree_.Branch(\"" << var->name_ << "\", &" << var->name_ << "_)),\n";
      }
      file << "  c_" << var->name_ << "_(false),\n";
    }
    file << "  " << com_end_2->name_ << "_(0),\n";
    file << "  b_" << com_end_2->name_ << "_(tree_.Branch(\"" << com_end_2->name_ << "\", &" << com_end_2->name_ << "_)),\n";
    file << "  c_" << com_end_2->name_ << "_(false){\n";
  }else{
    file << "  read_only_(false){\n";
  }
  file << "}\n\n";

  file << "small_tree::small_tree(const string &filename):\n";
  file << "  chain_(\"tree\",\"tree\"),\n";
  file << "  tree_(\"junk\",\"junk\"),\n";
  file << "  entry_(0),\n";
  if(com_vars.size()){
    const set<Variable>::const_iterator com_end_2 = --com_vars.end();
    file << "  read_only_(true),\n";
    for(set<Variable>::const_iterator var = com_vars.begin();
        var != com_end_2;
        ++var){
      file << "  " << var->name_ << "_(0),\n";
      if(Contains(var->type_, "vector")){
        file << "  p_" << var->name_ << "_(&" << var->name_ << "_),\n";
      }
      file << "  b_" << var->name_ << "_(NULL),\n";
      file << "  c_" << var->name_ << "_(false),\n";
    }
    file << "  " << com_end_2->name_ << "_(0),\n";
    if(Contains(com_end_2->type_, "vector")){
      file << "  p_" << com_end_2->name_ << "_(&" << com_end_2->name_ << "_),\n";
    }
    file << "  b_" << com_end_2->name_ << "_(NULL),\n";
    file << "  c_" << com_end_2->name_ << "_(false){\n";
  }else{
    file << "  read_only_(true){\n";
  }
  file << "  chain_.Add(filename.c_str());\n";
  for(set<Variable>::const_iterator var = com_vars.begin(); var != com_vars.end(); ++var){
    if(Contains(var->type_, "vector")){
      file << "  chain_.SetBranchAddress(\"" << var->name_ << "\", &p_" << var->name_ << "_, &b_" << var->name_ << "_);\n";
    }else{
      file << "  chain_.SetBranchAddress(\"" << var->name_ << "\", &" << var->name_ << "_, &b_" << var->name_ << "_);\n";
    }
  }
  file << "}\n\n";

  file << "void small_tree::Fill(){\n";
  file << "  if(read_only_){\n";
  file << "    throw std::logic_error(\"Trying to write to read-only tree\");\n";
  file << "  }else{\n";
  file << "    tree_.Fill();\n";
  file << "  }\n\n";

  file << "  //Resetting variables\n";
  for(set<Variable>::const_iterator var = com_vars.begin(); var != com_vars.end(); ++var){
    if(Contains(var->type_, "vector")){
      file << "  " << var->name_ << "_.clear();\n";
    }else if(Contains(var->type_, "tring")){
      file << "  " << var->name_ << "_ = \"\";\n";
    }else{
      file << "  " << var->name_ << "_ = static_cast<" << var->type_ << ">(bad_val_);\n";
    }
  }
  file << "}\n\n";

  file << "void small_tree::Write(){\n";
  file << "  if(read_only_){\n";
  file << "    throw std::logic_error(\"Trying to write to read-only tree.\");\n";
  file << "  }else{\n";
  file << "    tree_.Write();\n";
  file << "  }\n";
  file << "}\n\n";

  file << "string small_tree::Type() const{\n";
  file << "  return \"\";\n";
  file << "}\n\n";

  file << "small_tree::~small_tree(){\n";
  file << "}\n\n";

  file << "int small_tree::Add(const std::string &filename){\n";
  file << "  if(!read_only_){\n";
  file << "    throw std::logic_error(\"Trying to add files to tree opened for writing.\");\n";
  file << "  }\n";
  file << "  return chain_.Add(filename.c_str());\n";
  file << "}\n\n";

  file << "long small_tree::GetEntries() const{\n";
  file << "  if(read_only_){\n";
  file << "    return chain_.GetEntries();\n";
  file << "  }else{\n";
  file << "    return tree_.GetEntries();\n";
  file << "  }\n";
  file << "}\n\n";

  file << "void small_tree::GetEntry(const long entry){\n";
  file << "  if(!read_only_){\n";
  file << "    throw std::logic_error(\"Trying to read from write-only tree.\");\n";
  file << "  }\n\n";

  for(set<Variable>::const_iterator var = com_vars.begin(); var!= com_vars.end(); ++var){
    file << "  c_" << var->name_ << "_ = false;\n";
  }
  file << "  entry_ = chain_.LoadTree(entry);\n";
  file << "}\n\n";

  for(set<Variable>::const_iterator var = com_vars.begin(); var != com_vars.end(); ++var){
    file << var->type_ << " const & small_tree::" << var->name_ << "() const{\n";
    file << "  if(!read_only_){\n";
    file << "    throw std::logic_error(\"Trying to write to const tree.\");\n";
    file << "  }\n";
    file << "  if(!c_" << var->name_ << "_ && b_" << var->name_ <<"_){\n";
    file << "    b_" << var->name_ << "_->GetEntry(entry_);\n";
    file << "    c_" << var->name_ << "_ = true;\n";
    file << "  }\n";
    file << "  return " << var->name_ << "_;\n";
    file << "}\n\n";
  }

  for(set<Variable>::const_iterator var = com_vars.begin(); var != com_vars.end(); ++var){
    file << var->type_ << " & small_tree::" << var->name_ << "(){\n";
    file << "  if(read_only_ && !c_" << var->name_ << "_ && b_" << var->name_ <<"_){\n";
    file << "    b_" << var->name_ << "_->GetEntry(entry_);\n";
    file << "    c_" << var->name_ << "_ = true;\n";
    file << "  }\n";
    file << "  return " << var->name_ << "_;\n";
    file << "}\n\n";
  }

  for(set<Variable>::const_iterator var = all_vars.begin(); var != all_vars.end(); ++var){
    if(com_vars.find(*var) != com_vars.end()) continue;
    file << var->type_ << " const & small_tree::" << var->name_ << "() const{\n";
    file << "  throw std::logic_error(\"" << var->name_
         << " does not exist in this small_tree version.\");\n";
    file << "}\n\n";
  }

  for(set<Variable>::const_iterator var = all_vars.begin(); var != all_vars.end(); ++var){
    if(com_vars.find(*var) != com_vars.end()) continue;
    file << var->type_ << " & small_tree::" << var->name_ << "(){\n";
    file << "  throw std::logic_error(\"" << var->name_
         << " does not exist in this small_tree version.\");\n";
    file << "}\n\n";
  }

  for(size_t i = 0; i < names.size(); ++i){
    file << "#include \"small_tree_" << names.at(i) << ".hpp\"\n";
  }
  file << "small_tree* NewTree(const std::type_info &type){\n\n";
  file << "  if(type == typeid(small_tree)) return new small_tree;\n";
  for(size_t i = 0; i < names.size(); ++i){
    file << "  else if(type == typeid(small_tree_" << names.at(i) << ")) return static_cast<small_tree*>(new small_tree_" << names.at(i) << ");\n";
  }
  file << "  else return new small_tree;\n";
  file << "}\n\n";

  file.close();
}

void WriteSepHeader(const pair<string, set<Variable> > &sep_vars){
  string name = sep_vars.first;
  string NAME = ToCaps(name);
  set<Variable> vars = sep_vars.second;
  ofstream file(("inc/small_tree_"+name+".hpp").c_str());

  file << "// small_tree_" << name << ": " << name << " version of small_tree to handle reduce tree ntuples\n";
  file << "// File generated with generate_small_tree.exe\n\n";

  file << "#ifndef H_SMALL_TREE_" << NAME << "\n";
  file << "#define H_SMALL_TREE_" << NAME << "\n\n";

  file << "#include <vector>\n";
  file << "#include <string>\n\n";

  file << "#include \"TTree.h\"\n";
  file << "#include \"TChain.h\"\n\n";

  file << "#include \"small_tree.hpp\"\n\n";

  file << "class small_tree_" << name << " : public small_tree{\n";
  file << "public:\n";
  file << "  small_tree_" << name << "(); // Constructor to create tree\n";
  file << "  small_tree_" << name << "(const std::string &filename); // Constructor to read tree\n\n";

  file << "  virtual void GetEntry(const long entry);\n\n";

  file << "  virtual void Fill();\n\n";

  file << "  virtual std::string Type() const;\n\n";

  file << "  virtual ~small_tree_" << name << "();\n\n";

  for(set<Variable>::const_iterator var = vars.begin();
      var != vars.end();
      ++var){
    file << "  virtual " << var->type_ << " const & " << var->name_ << "() const;\n";
    file << "  virtual " << var->type_ << " & " << var->name_ << "();\n";
  }
  file << '\n';

  file << "private:\n";
  for(set<Variable>::const_iterator var = vars.begin();
      var != vars.end();
      ++var){
    file << "  " << var->type_ << ' ' << var->name_ << "_;\n";
    if(Contains(var->type_, "vector")){
      file << "  " << var->type_ << " *p_" << var->name_ << "_;\n";
    }
    file << "  TBranch *b_" << var->name_ << "_;\n";
    file << "  mutable bool c_" << var->name_ << "_;\n";
  }
  file << "};\n\n";

  file << "#endif" << endl;

  file.close();
}

void WriteSepSource(const pair<string, set<Variable> > &sep_vars){
  string name = sep_vars.first;
  string NAME = ToCaps(name);
  set<Variable> vars = sep_vars.second;
  ofstream file(("src/small_tree_"+name+".cpp").c_str());

  file << "// small_tree_" << name << ": " << name << " version of small_tree to handle reduce tree ntuples\n";
  file << "//File generated with generate_small_tree.exe\n\n";

  file << "#include \"small_tree.hpp\"\n\n";
  file << "#include \"small_tree_" << name << ".hpp\"\n\n";

  file << "#include <stdexcept>\n";
  file << "#include <string>\n";
  file << "#include <vector>\n\n";

  file << "#include \"TTree.h\"\n";
  file << "#include \"TBranch.h\"\n";
  file << "#include \"TChain.h\"\n\n";

  file << "using namespace std;\n\n";

  file << "small_tree_" << name << "::small_tree_" << name << "():\n";
  if(vars.size()){
    const set<Variable>::const_iterator vars_end_2 = --vars.end();
    file << "  small_tree(),\n";
    for(set<Variable>::const_iterator var = vars.begin();
        var != vars_end_2;
        ++var){
      file << "  " << var->name_ << "_(0),\n";
      if(Contains(var->type_, "vector")){
        file << "  p_" << var->name_ << "_(&" << var->name_ << "_),\n";
        file << "  b_" << var->name_ << "_(tree_.Branch(\"" << var->name_ << "\", &p_" << var->name_ << "_)),\n";
      }else{
        file << "  b_" << var->name_ << "_(tree_.Branch(\"" << var->name_ << "\", &" << var->name_ << "_)),\n";
      }
      file << "  c_" << var->name_ << "_(false),\n";
    }
    file << "  " << vars_end_2->name_ << "_(0),\n";
    if(Contains(vars_end_2->type_, "vector")){
      file << "  p_" << vars_end_2->name_ << "_(&" << vars_end_2->name_ << "_),\n";
      file << "  b_" << vars_end_2->name_ << "_(tree_.Branch(\"" << vars_end_2->name_ << "\", &p_" << vars_end_2->name_ << "_)),\n";
    }else{
      file << "  b_" << vars_end_2->name_ << "_(tree_.Branch(\"" << vars_end_2->name_ << "\", &" << vars_end_2->name_ << "_)),\n";
    }
    file << "  c_" << vars_end_2->name_ << "_(false){\n";
  }else{
    file << "  small_tree(){\n";
  }
  file << "}\n\n";

  file << "small_tree_" << name << "::small_tree_" << name << "(const string &filename):\n";
  if(vars.size()){
    const set<Variable>::const_iterator vars_end_2 = --vars.end();
    file << "  small_tree(filename),\n";
    for(set<Variable>::const_iterator var = vars.begin();
        var != vars_end_2;
        ++var){
      file << "  " << var->name_ << "_(0),\n";
      if(Contains(var->type_, "vector")){
        file << "  p_" << var->name_ << "_(&" << var->name_ << "_),\n";
      }
      file << "  b_" << var->name_ << "_(NULL),\n";
      file << "  c_" << var->name_ << "_(false),\n";
    }
    file << "  " << vars_end_2->name_ << "_(0),\n";
    if(Contains(vars_end_2->type_, "vector")){
      file << "  p_" << vars_end_2->name_ << "_(&" << vars_end_2->name_ << "_),\n";
    }
    file << "  b_" << vars_end_2->name_ << "_(NULL),\n";
    file << "  c_" << vars_end_2->name_ << "_(false){\n";
  }else{
    file << "  small_tree(filename){\n";
  }
  for(set<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    if(Contains(var->type_, "vector")){
      file << "  chain_.SetBranchAddress(\"" << var->name_ << "\", &p_" << var->name_ << "_, &b_" << var->name_ << "_);\n";
    }else{
      file << "  chain_.SetBranchAddress(\"" << var->name_ << "\", &" << var->name_ << "_, &b_" << var->name_ << "_);\n";
    }
  }
  file << "}\n\n";

  file << "void small_tree_" << name << "::Fill(){\n";
  file << "  small_tree::Fill();\n";

  file << "  //Resetting variables\n";
  for(set<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    if(Contains(var->type_, "vector")){
      file << "  " << var->name_ << "_.clear();\n";
    }else if(Contains(var->type_, "tring")){
      file << "  " << var->name_ << "_ = \"\";\n";
    }else{
      file << "  " << var->name_ << "_ = static_cast<" << var->type_ << ">(bad_val_);\n";
    }
  }
  file << "}\n\n";

  file << "string small_tree_" << name << "::Type() const{\n";
  file << "  return \"" << name << "\";\n";
  file << "}\n\n";

  file << "small_tree_" << name << "::~small_tree_" << name << "(){\n";
  file << "}\n\n";

  file << "void small_tree_" << name << "::GetEntry(const long entry){\n";
  file << "  small_tree::GetEntry(entry);\n\n";

  for(set<Variable>::const_iterator var = vars.begin(); var!= vars.end(); ++var){
    file << "  c_" << var->name_ << "_ = false;\n";
  }
  file << "}\n\n";

  for(set<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    file << var->type_ << " const & small_tree_" << name << "::" << var->name_ << "() const{\n";
    file << "  if(!read_only_){\n";
    file << "    throw std::logic_error(\"Trying to write to const tree.\");\n";
    file << "  }\n";
    file << "  if(!c_" << var->name_ << "_ && b_" << var->name_ <<"_){\n";
    file << "    b_" << var->name_ << "_->GetEntry(entry_);\n";
    file << "    c_" << var->name_ << "_ = true;\n";
    file << "  }\n";
    file << "  return " << var->name_ << "_;\n";
    file << "}\n\n";
  }

  for(set<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    file << var->type_ << " & small_tree_" << name << "::" << var->name_ << "(){\n";
    file << "  if(read_only_ && !c_" << var->name_ << "_ && b_" << var->name_ <<"_){\n";
    file << "    b_" << var->name_ << "_->GetEntry(entry_);\n";
    file << "    c_" << var->name_ << "_ = true;\n";
    file << "  }\n";
    file << "  return " << var->name_ << "_;\n";
    file << "}\n\n";
  }

  file.close();
}

void GenerateEventHandlerBaseHeader(){
  ofstream file("inc/event_handler_base.hpp");

  file << " // event_handler_base: base class for reduced tree production\n";
  file << " // File generated with generate_small_tree.exe\n\n";

  file << "#ifndef H_EVENT_HANDLER_BASE\n";
  file << "#define H_EVENT_HANDLER_BASE\n\n";

  file << "#include <string>\n\n";

  file << "#include \"TString.h\"\n\n";

  file << "#include \"phys_objects.hpp\"\n\n";

  file << "class event_handler_base : public phys_objects{\n";
  file << "public:\n";
  file << "  event_handler_base(const std::string &file_name);\n\n";

  file << "  virtual void ReduceTree(int num_entries,\n";
  file << "                          const TString &out_file_name,\n";
  file << "                          int num_total_entries) = 0;\n\n";

  file << "  virtual ~event_handler_base();\n";

  file << "};\n\n";

  file << "#endif\n";

  file.close();
}

void GenerateEventHandlerBaseSource(){
  ofstream file("src/event_handler_base.cpp");

  file << "// event_handle_base: base class for reduced tree production\n";
  file << "// File generated with generate_small_tree.exe\n\n";

  file << "#include \"event_handler_base.hpp\"\n\n";

  file << "#include <string>\n";

  file << "using namespace std;\n\n";

  file << "event_handler_base::event_handler_base(const string &file_name):\n";
  file << "  phys_objects(file_name){\n";
  file << "}\n\n";

  file << "event_handler_base::~event_handler_base(){\n";
  file << "}\n\n";
}

void GenerateEventHandlerHeader(const string &name){
  string NAME = ToCaps(name);
  ofstream file(("inc/event_handler_"+name+".hpp").c_str());

  file << " // event_handler_" << name << ": derived class for specialized reduced tree production\n";
  file << " // File generated with generate_small_tree.exe\n\n";

  file << "#ifndef H_EVENT_HANDLER_" << NAME << "\n";
  file << "#define H_EVENT_HANDLER_" << NAME << "\n";

  file << "#include <string>\n";

  file << "#include \"TString.h\"\n\n";

  file << "#include \"event_handler_base.hpp\"\n\n";

  file << "class event_handler_" << name << " : public event_handler_base{\n";
  file << "public:\n";
  file << "  event_handler_" << name << "(const std::string &file_name);\n\n";

  file << "  virtual void ReduceTree(int num_entries,\n";
  file << "                          const TString &out_file_name,\n";
  file << "                          int num_total_entries);\n\n";

  file << "  virtual ~event_handler_" << name << "();\n";

  file << "};\n\n";

  file << "#endif\n";

  file.close();
}

void GenerateEventHandlerSource(const string &name){
  ofstream file(("src/event_handler_"+name+".cpp").c_str());

  file << " // event_handler_" << name << ": derived class for specialized reduced tree production\n";
  file << " // File generated with generate_small_tree.exe\n\n";

  file << "#include \"event_handler_" << name << ".hpp\"\n\n";

  file << "#include <string>\n";

  file << "#include \"TString.h\"\n\n";

  file << "#include \"event_handler_base.hpp\"\n\n";

  file << "using namespace std;\n\n";

  file << "event_handler_" << name << "::event_handler_" << name << "(const string &file_name):\n";
  file << "  event_handler_base(file_name){\n";
  file << "}\n\n";

  file << "void event_handler_" << name << "::ReduceTree(int /*num_entries*/, const TString &/*out_file_name*/, int /*num_total_entries*/){\n";
  file << "}\n\n";

  file << "event_handler_" << name << "::~event_handler_" << name << "(){\n";
  file << "}\n";

  file.close();
}

void GenerateEventHandlerMergeHeader(const vector<string> &names){
  ofstream file("inc/event_handler.hpp");

  file << " // event_handler: main class for reduced tree production\n";
  file << " //File generated with generate_small_tree.exe\n\n";

  file << "#ifndef H_EVENT_HANDLER\n";
  file << "#define H_EVENT_HANDLER\n\n";

  file << "#include <string>\n";

  file << "#include \"TString.h\"\n\n";

  file << "#include \"event_handler_base.hpp\"\n\n";
  for(vector<string>::const_iterator name = names.begin(); name != names.end(); ++name){
    file << "#include \"event_handler_" << *name << ".hpp\"\n";
  }

  file << "class event_handler{\n";
  file << "public:\n";
  file << "  event_handler(const std::string &file_name, const std::string &type);\n";
  file << "  void ReduceTree(int num_entries, const TString &out_file_name, int num_total_entries);\n";
  file << "  long TotalEntries() const;\n";
  file << "  short GetVersion() const;\n";
  file << "  const std::string& SampleName() const;\n";
  file << "  const std::string& SampleName(const std::string &sample_name);\n";
  file << "  void SetFile(const std::string &file, bool is_8TeV = false);\n";
  file << "  void AddFiles(const std::string &file);\n";
  file << "  ~event_handler();\n";
  file << "private:\n";
  file << "  event_handler_base *ehb;\n";
  file << "  static event_handler_base * LookUpType(const std::string &file_name, const std::string &type);\n";
  file << "};\n\n";

  file << "#endif\n";

  file.close();
}

void GenerateEventHandlerMergeSource(const vector<string> &names){
  ofstream file("src/event_handler.cpp");

  file << "// event_handler: main class for reduced tree production\n";
  file << " //File generated with generate_small_tree.exe\n\n";

  file << "#include \"event_handler.hpp\"\n\n";

  file << "#include <typeinfo>\n";
  file << "#include <string>\n\n";

  file << "#include \"TString.h\"\n\n";

  file << "#include \"event_handler_base.hpp\"\n\n";
  for(vector<string>::const_iterator name = names.begin(); name != names.end(); ++name){
    file << "#include \"event_handler_" << *name << ".hpp\"\n";
  }

  file << "using namespace std;\n\n";

  file << "#include \"event_handler_base.hpp\"\n\n";
  for(vector<string>::const_iterator name = names.begin(); name != names.end(); ++name){
    file << "#include \"event_handler_" << *name << ".hpp\"\n";
  }

  file << "event_handler::event_handler(const string &file_name, const string &type):\n";
  file << "  ehb(LookUpType(file_name, type)){\n";
  file << "}\n\n";

  file << "void event_handler::ReduceTree(int num_entries, const TString &out_file_name, int num_total_entries){\n";
  file << "  ehb->ReduceTree(num_entries, out_file_name, num_total_entries);\n";
  file << "}\n\n";

  if(names.size()){
    file << "event_handler_base * event_handler::LookUpType(const string &file_name, const string &type){\n";
    file << "  if(type ==  \"" << names.front() << "\"){\n";
    file << "    return new event_handler_" << names.front() << "(file_name);\n";
    for(size_t itype = 1; itype < names.size(); ++itype){
      file << "  }else if(type == \"" << names.at(itype) << "\"){\n";
      file << "    return new event_handler_" << names.at(itype) << "(file_name);\n";
    }
    file << "  }else{\n";
    file << "    return NULL;\n";
    file << "  }\n";
  }else{
    file << "event_handler_base * LookUpType(const string &/*file_name*/, const string &/*type*/){\n";
    file << "  return NULL;\n";
  }
  file << "}\n\n";

  file << "long event_handler::TotalEntries() const{\n";
  file << "  return ehb->TotalEntries();\n";
  file << "}\n\n";

  file << "short event_handler::GetVersion() const{\n";
  file << "  return ehb->GetVersion();\n";
  file << "}\n\n";

  file << "const std::string& event_handler::SampleName() const{\n";
  file << "  return ehb->SampleName();\n";
  file << "}\n\n";

  file << "const std::string& event_handler::SampleName(const std::string &sample_name){\n";
  file << "  return ehb->SampleName(sample_name);\n";
  file << "}\n\n";

  file << "void event_handler::SetFile(const std::string &file, bool is_8TeV){\n";
  file << "  ehb->SetFile(file, is_8TeV);\n";
  file << "}\n\n";

  file << "void event_handler::AddFiles(const std::string &file){\n";
  file << "  ehb->AddFiles(file);\n";
  file << "}\n\n";

  file << "event_handler::~event_handler(){\n";
  file << "  if(ehb){\n";
  file << "    delete ehb;\n";
  file << "    ehb = NULL;\n";
  file << "  }\n";
  file << "}\n\n";
}
