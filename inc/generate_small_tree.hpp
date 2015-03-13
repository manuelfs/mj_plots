#ifndef H_GENERATE_SMALL_TREE
#define H_GENERATE_SMALL_TREE

#include <vector>
#include <set>
#include <string>

class Variable{
public:
  Variable():
    type_(""),
    name_(""){
  }

  Variable(const std::string &type,
           const std::string &name):
    type_(type),
    name_(name){
  }

  bool operator<(const Variable& var) const{
    return type_<var.type_ || (type_==var.type_ && name_<var.name_);
  }

  std::string type_, name_;
};

bool Contains(const std::string &text, const std::string &pattern);

void WriteBaseHeader(const std::set<Variable> &all_vars,
                     const std::set<Variable> &com_vars,
                     const std::vector<std::string> &names);
void WriteBaseSource(const std::set<Variable> &all_vars,
                     const std::set<Variable> &com_vars,
                     const std::vector<std::string> &names);
void WriteSepHeader(const std::pair<std::string, std::set<Variable> > &vars);
void WriteSepSource(const std::pair<std::string, std::set<Variable> > &vars);

void GenerateEventHandlerBaseHeader();
void GenerateEventHandlerBaseSource();
void GenerateEventHandlerHeader(const std::string& name);
void GenerateEventHandlerSource(const std::string& name);
void GenerateEventHandlerMergeHeader(const std::vector<std::string> &names);
void GenerateEventHandlerMergeSource(const std::vector<std::string> &names);

#endif
