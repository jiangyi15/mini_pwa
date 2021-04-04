
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnUserParameters.h"

#include "Minuit2/FCNBase.h"

#include "decay.h"

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

class AmplitudeModel : public FCNBase {
public:
  AmplitudeModel(DecayGroup *decay_group)
      : decay_group(decay_group), vm(), data(), phsp() {
    this->decay_group->init_params(&(this->vm));
  };
  DecayGroup *decay_group;
  VarManager vm;
  std::map<std::string, ChainData *> data;
  std::map<std::string, ChainData *> phsp;
  size_t n_data, n_phsp;
  void set_data(size_t n, std::map<std::string, ChainData *> &data) {
    this->n_data = n;
    this->data = data;
  };
  void set_phsp(size_t n, std::map<std::string, ChainData *> &phsp) {
    this->n_phsp = n;
    this->phsp = phsp;
  };

  double operator()(const std::vector<double> &par) const {
    this->vm.set_params(par);
    auto amp = this->decay_group->get_amp2s(this->n_data, this->data);
    auto amp_mc = this->decay_group->get_amp2s(this->n_phsp, this->phsp);
    double nll = 0.0;
    for (int i = 0; i < this->n_data; i++) {
      nll += -log(amp.ptr[i]);
    }
    double int_mc = 0.0;
    for (int i = 0; i < this->n_phsp; i++) {
      int_mc += amp_mc.ptr[i];
    }
    return nll + this->n_data * log(int_mc / this->n_phsp);
  }

  FunctionMinimum minimize() {
    MnUserParameters upar;
    for (auto i : this->vm.vars) {
      upar.Add(i.first, i.second, 0.1);
    }
    MnMigrad migrad(*this, upar);
    FunctionMinimum min = migrad();
    return min;
  }
  virtual double Up() const override { return 1.0; };
};
