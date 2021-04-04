#include "data.h"
#include "decay.h"
#include "variable.h"
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "fcn.h"

int main(int argc, char **argv) {

  auto a = new Particle("a");
  auto b = new Particle("b");
  auto c = new Particle("c");

  auto x = new Tensor<double>({100});
  for (int i = 0; i < 100; i++) {
    x->ptr[i] = cos(i) * cos(i);
  }

  auto y = new Tensor<double>({100});
  for (int i = 0; i < 100; i++) {
    y->ptr[i] = sin(i) * sin(i);
  }
  auto d2 = new DecayData({x, x, x});
  auto d4 = new ChainData({{d2}, {{"a", x}}});

  auto data = std::map<std::string, ChainData *>({
      {"[a->b+c]", d4},
  });
  auto d2_p = new DecayData({y, y, y});
  auto d4_p = new ChainData({{d2_p}, {{"a", y}}});

  auto phsp = std::map<std::string, ChainData *>({
      {"[a->b+c]", d4_p},
  });

  auto dec = new Decay(a, {b, c});
  auto decs = new DecayChain({
      dec,
  });

  auto dg = new DecayGroup({
      decs,
  });

  AmplitudeModel amp(dg);

  amp.vm.set_params({{"a_mass", 2.0}, {"a_width", 1.0}});

  amp.set_data(100, data);
  amp.set_phsp(100, phsp);

  auto min = amp.minimize();
  std::cout << min << std::endl;
  delete a, b, c, x;
  delete d2, d4, dec, decs, dg;
  return 0;
}
