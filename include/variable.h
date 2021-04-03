#include <functional>
#include <map>
#include <string>
#include <fstream>
#ifndef VARIABLE_H_
#define VARIABLE_H_ 1

void save_params(std::string file_name, std::map<std::string, double> params) {
  std::ofstream out(file_name);
  out << "{" << std::endl;
  for (auto iter : params) {
    out << "\t\"" << iter.first << "\": \"" << iter.second << "\","
        << std::endl;
  }
  out << "}";
  out.close();
}

class VarManager {
public:
  VarManager() : vars(){};
  std::map<std::string, double> vars;
  std::function<double()> add_var(std::string name) {
    auto iter = this->vars.find("123");
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

  void save_params(std::string file_name) {
    ::save_params(file_name, this->vars);
  }
};

#endif
