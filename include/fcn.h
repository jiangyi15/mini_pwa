
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnUserParameters.h"

#include "Minuit2/FCNBase.h"

#ifndef M_PWA_FCN_H
#define M_PWA_FCN_H

#include "decay.h"

using namespace ROOT::Minuit2;

class AmplitudeModel : public FCNBase {
public:
  AmplitudeModel(DecayGroup *decay_group)
      : decay_group(decay_group), vm(), fix_params(), data(), phsp() {
    this->decay_group->init_params(&(this->vm));
  };
  DecayGroup *decay_group;
  mutable VarManager vm;
  std::map<std::string, double> fix_params;
  std::map<std::string, ChainData *> data;
  Tensor<double> *data_weight = nullptr;
  std::map<std::string, ChainData *> phsp;
  Tensor<double> *phsp_weight = nullptr;
  size_t n_data, n_phsp;
  void set_data(size_t n, std::map<std::string, ChainData *> &data) {
    this->n_data = n;
    this->data = data;
  };
  void set_phsp(size_t n, std::map<std::string, ChainData *> &phsp) {
    this->n_phsp = n;
    this->phsp = phsp;
  };

  double operator()(const std::vector<double> &par) const override;

  double get_n_phsp() const;

  double get_integral() const;

  FunctionMinimum minimize();
  virtual double Up() const override { return 1.0; };
};

#endif
