
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnUserParameters.h"

#include "Minuit2/FCNBase.h"

#include "decay.h"

using namespace ROOT::Minuit2;

class AmplitudeModel : public FCNBase {
public:
  AmplitudeModel(DecayGroup *decay_group)
      : decay_group(decay_group), vm(), fix_params(), data(), phsp() {
    this->decay_group->init_params(&(this->vm));
  };
  DecayGroup *decay_group;
  VarManager vm;
  std::map<std::string, double> fix_params;
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

  double operator()(const std::vector<double> &par) const override {
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
      if (this->fix_params.find(i.first) != this->fix_params.end()) {
      } else {
        upar.Add(i.first, i.second, 0.1);
      }
    }
    MnMigrad migrad(*this, upar);
    FunctionMinimum min = migrad();
    return min;
  }
  virtual double Up() const override { return 1.0; };
};
