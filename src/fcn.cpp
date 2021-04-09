#include "fcn.h"

double AmplitudeModel::operator()(const std::vector<double> &par) const {
  this->vm.set_params(par);
  auto amp = this->decay_group->get_amp2s(this->n_data, this->data);
  double nll = 0.0;
  if (this->data_weight == nullptr) {
    for (int i = 0; i < this->n_data; i++) {
      nll += -log(amp.ptr[i]);
    }
  } else {
    for (int i = 0; i < this->n_data; i++) {
      nll += -data_weight->ptr[i] * log(amp.ptr[i]);
    }
  }
  double int_mc = this->get_integral();
  return nll + this->n_data * log(int_mc / this->get_n_phsp());
}

double AmplitudeModel::get_n_phsp() const {
  if (this->phsp_weight == nullptr) {
    return this->n_phsp;
  } else {
    double sw = 0.0;
    for (int i = 0; i < this->n_phsp; i++) {
      sw += this->phsp_weight->ptr[i];
    }
    return sw;
  }
}

double AmplitudeModel::get_integral() const {
  auto amp_mc = this->decay_group->get_amp2s(this->n_phsp, this->phsp);
  double int_mc = 0.0;
  if (this->phsp_weight == nullptr) {
    for (int i = 0; i < this->n_phsp; i++) {
      int_mc += amp_mc.ptr[i];
    }
  } else {
    for (int i = 0; i < this->n_phsp; i++) {
      int_mc += this->phsp_weight->ptr[i] * amp_mc.ptr[i];
    }
  }
  return int_mc;
}

FunctionMinimum AmplitudeModel::minimize() {
  MnUserParameters upar;
  for (auto i : this->vm.vars) {
    if (this->fix_params.find(i.first) != this->fix_params.end()) {
      upar.Add(i.first, this->fix_params[i.first]);
    } else {
      if (this->bound.find(i.first) != this->bound.end()) {
        auto bound = this->bound[i.first];
        upar.Add(i.first, i.second, 0.1, bound.first, bound.second);
      } else {
        upar.Add(i.first, i.second, 0.1);
      }
    }
  }
  MnMigrad migrad(*this, upar);
  FunctionMinimum min = migrad();
  return min;
}
