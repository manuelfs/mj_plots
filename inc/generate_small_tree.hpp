#ifndef H_GENERATE_SMALL_TREE
#define H_GENERATE_SMALL_TREE

#include <vector>
#include <string>

struct Variable{
  std::string name, type;
};

std::vector<Variable> ParseVariables(const std::vector<std::string> &variables);
bool Contains(const std::string &text, const std::string &pattern);

#endif
