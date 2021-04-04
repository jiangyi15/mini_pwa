#include "data.h"
#include "decay.h"
#include "variable.h"
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnUserParameters.h"

#include "Minuit2/FCNBase.h"

using namespace ROOT::Minuit2;

class MyFCN : public FCNBase {
public:
  double operator()(const std::vector<double> &par) const {

    double x = par[0];
    double y = par[1];
    double z = par[2];
    double w = par[3];
    double x0 = par[4];
    double y0 = par[5];
    double z0 = par[6];
    double w0 = par[7];

    return ((1. / 70.) * (21 * x * x + 20 * y * y + 19 * z * z - 14 * x * z -
                          20 * y * z) +
            w * w +
            (1. / 70.) * (21 * x0 * x0 + 20 * y0 * y0 + 19 * z0 * z0 -
                          14 * x0 * z0 - 20 * y0 * z0) +
            w0 * w0);
  }

  double Up() const { return 1.; }
};

void minimize() {
  auto fcn = MyFCN();
  MnUserParameters upar;
  upar.Add("x", 1., 0.1);
  upar.Add("y", 1., 0.1);
  upar.Add("z", 1., 0.1);
  upar.Add("w", 1., 0.1);
  upar.Add("x0", 1., 0.1);
  upar.Add("y0", 1., 0.1);
  upar.Add("z0", 1., 0.1);
  upar.Add("w0", 1., 0.1);

  MnMigrad migrad(fcn, upar);
  FunctionMinimum min = migrad();
  std::cout << min << std::endl;
}

int main(int argc, char **argv) {
  minimize();
  std::map<std::string, double> params;
  params["s"] = 2;
  save_params("a.json", params);

  auto a = new Particle("a");
  auto b = new Particle("b");
  auto c = new Particle("c");

  auto x = new Tensor<double>({100});
  for (int i = 0; i < 100; i++) {
    x->ptr[i] = 1;
  }
  auto d2 = new DecayData({x, x, x});
  auto d4 = new ChainData({{d2}, {{"a", x}}});

  auto dec = new Decay(a, {b, c});
  auto decs = new DecayChain({
      dec,
  });

  auto dg = new DecayGroup({
      decs,
  });

  auto vm = VarManager();
  dg->init_params(&vm);

  vm.set_params({{"a_mass", 2.0}, {"a_width", 1.0}});

  std::cout << dg->get_amp2s(1,
                             {
                                 {"[a->b+c]", d4},
                             })
                   .ptr[0]
            << std::endl;

  vm.save_params("b.json");
  return 0;
}
