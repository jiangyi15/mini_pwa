#include "variable.h"

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
