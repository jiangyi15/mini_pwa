#include "data.h"
#include "decay.h"
#include "variable.h"
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "dfunction.h"
#include "fcn.h"

#include "cg_coeff.h"

int main(int argc, char **argv) {

  std::cout << "cg_coeff:" << cg_coeff2(2, 1, -2, 1, 1, -1) << std::endl;

  std::cout << "start";
  auto w = small_d_weight(2);
  std::cout << "end";

  auto a = new Particle("a");
  auto b = new Particle("b");
  auto c = new Particle("c");
  auto d = new Particle("d");
  auto r = new Particle("r", 1, -1);

  auto dec1 = new Decay(a, {r, d});
  auto dec2 = new Decay(r, {b, c});
  auto decs = new DecayChain({dec1, dec2});

  auto dg = new DecayGroup({
      decs,
  });

  auto x = new Tensor<double>({100});
  for (int i = 0; i < x->total_shape(); i++) {
    x->ptr[i] = sin(i);
  }

  auto dfun = small_d_function(2, &w, x);

  auto y = new Tensor<double>({100});
  for (int i = 0; i < 1; i++) {
    y->ptr[i] = 1.0;
  }
  auto d1 = new EularAngle(x, x, x);
  auto d2 = new DecayData(d1, x);
  auto d4 = new ChainData({d2, d2}, {{"a", x}, {"r", x}});

  auto data = std::map<std::string, ChainData *>({
      {"[a->r+d,r->b+c]", d4},
  });
  auto d1_p = new EularAngle(y, y, y);

  auto d2_p = new DecayData(d1_p, y);
  auto d4_p = new ChainData({d2_p, d2_p}, {{"a", y}, {"r", y}});

  auto phsp = std::map<std::string, ChainData *>({
      {"[a->r+d,r->b+c]", d4_p},
  });

  AmplitudeModel amp(dg);

  amp.vm.set_params({{"r_mass", 2.0}, {"r_width", 1.0}});
  amp.fix_params["r_width"] = 1.0;

  amp.set_data(x->shape[0], data);
  amp.set_phsp(y->shape[0], phsp);

  auto min = amp.minimize();
  std::cout << min << std::endl;
  delete a, b, c, d, r, x;
  delete d2, d4, dec1, dec2, decs, dg;
  return 0;
}
