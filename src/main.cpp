#include "data.h"
#include "decay.h"
#include "variable.h"
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  std::map<std::string, double> params;
  params["s"] = 2;
  save_params("a.json", params);

  auto a = Particle("a");
  auto b = Particle("b");
  auto c = Particle("c");

  auto x = new Tensor<double>({100});
  for (int i = 0; i < 100; i++) {
    x->ptr[i] = 1;
  }
  auto d2 = new DecayData({x, x, x});
  auto d4 = new ChainData({{d2}, {{"a", x}}});

  auto dec = Decay(&a, {&b, &c});
  auto decs = DecayChain({
      &dec,
  });
  auto dg = DecayGroup({
      &decs,
  });

  auto vm = VarManager();
  dg.init_params(&vm);

  vm.set_params({{"a_mass", 2.0}, {"a_width", 1.0}});

  std::cout << dg.get_amp2s({
                   {"[a->b+c]", d4},
               })
            << std::endl;

  vm.save_params("b.json");
  return 0;
}
