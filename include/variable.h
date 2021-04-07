#include <fstream>
#include <functional>
#include <map>
#include <string>
#include <vector>
#ifndef VARIABLE_H_
#define VARIABLE_H_ 1

void save_params(std::string file_name, std::map<std::string, double> params);

class VarManager {
public:
  VarManager() : vars(){};
  std::map<std::string, double> vars;
  std::function<double()> add_var(std::string name) {
    auto iter = this->vars.find(name);
    if (iter == vars.end()) {
      this->vars[name] = 1.0;
    }
    return [this, name]() { return this->vars.at(name); };
  }

  void set_params(std::map<std::string, double> new_vars) {
    for (auto i : new_vars) {
      auto iter = this->vars.find(i.first);
      if (iter != vars.end()) {
        this->vars[i.first] = i.second;
      }
    }
  }
  void set_params(std::vector<double> new_vars) {
    int idx = 0;
    for (auto i : this->vars) {
      this->vars[i.first] = new_vars[idx++];
    }
  }

  void save_params(std::string file_name) {
    ::save_params(file_name, this->vars);
  }
};

#endif
