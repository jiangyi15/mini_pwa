#include "data.h"
#include "variable.h"
#include <complex>
#include <functional>
#include <iostream>
#include <string>
#include <vector>
#ifndef M_PWADECAY_H_
#define M_PWADECAY_H_ 1

class BaseParticle {
public:
  int J, P;
  BaseParticle(std::string name, int J, int P) : name(name), J(J), P(P){};
  BaseParticle(std::string name) : name(name), J(0), P(-1){};
  std::string name;
  virtual void init_params(VarManager *vm) {}
  Tensor<std::complex<double>> get_amp(size_t n, Tensor<double> *data) {
    Tensor<std::complex<double>> ret(n);
    for (int i = 0; i < n; i++) {
      ret.ptr[i] = 1.;
    }
    return ret;
  };
  std::string to_string() { return this->name; }
};

const std::function<double()> zeros_f = [] { return 0.0; };

class Particle : public BaseParticle {
public:
  Particle(std::string name) : BaseParticle(name){};
  std::function<double()> mass = zeros_f;
  std::function<double()> width = zeros_f;
  virtual void init_params(VarManager *vm) {
    this->mass = vm->add_var(this->to_string() + "_mass");
    this->width = vm->add_var(this->to_string() + "_width");
  }

  Tensor<std::complex<double>> get_amp(size_t n, Tensor<double> *data) {
    Tensor<std::complex<double>> ret(n);
    for (int i = 0; i < n; i++) {
      double m = data->ptr[i];
      double m0 = this->mass();
      double g0 = this->width();
      ret.ptr[i] = 1. / std::complex(m0 * m0 - m * m, -m0 * g0);
    }
    return ret;
  };
};

class BaseDecay {
public:
  BaseDecay(BaseParticle *core, std::vector<BaseParticle *> decay)
      : core(core), decay(decay){};
  BaseParticle *core;
  std::vector<BaseParticle *> decay;
  virtual void init_params(VarManager *vm) {}
  virtual Tensor<std::complex<double>> get_amp(size_t n, DecayData *data) {
    auto ret = Tensor<std::complex<double>>(n);
    for (int i = 0; i < n; i++) {
      ret.ptr[i] = data->angle->alpha->ptr[i];
    }
    return ret;
  };
  std::string to_string() {
    std::ostringstream ostr;
    ostr << this->core->to_string() << "->";
    int idx = 0;
    for (auto i : this->decay) {
      if (idx == 0) {
        ostr << i->to_string();
      } else {
        ostr << "+" << i->to_string();
      }
      idx++;
    }
    return ostr.str();
  }
};

class Decay : public BaseDecay {
public:
  Decay(BaseParticle *core, std::vector<BaseParticle *> decay)
      : BaseDecay(core, decay){};
};

class BaseDecayChain {
public:
  BaseDecayChain(std::vector<BaseDecay *> decs) : decays(decs) {
    this->top = decs[0]->core;
  };
  BaseParticle *top;
  std::vector<BaseDecay *> decays;
  std::function<double()> total_r = zeros_f;
  std::function<double()> total_i = zeros_f;
  virtual void init_params(VarManager *vm) {
    this->total_r = vm->add_var(this->to_string() + "total_r");
    this->total_i = vm->add_var(this->to_string() + "total_i");
    for (auto i : this->decays) {
      i->init_params(vm);
      if (i->core != this->top) {
        i->core->init_params(vm);
      }
    }
  }

  virtual Tensor<std::complex<double>> get_amp(size_t n, ChainData *data) {
    Tensor<std::complex<double>> ret(n);
    for (int i = 0; i < n; i++) {
      ret.ptr[i] = std::complex(this->total_r(), this->total_i());
    }
    for (int idx = 0; idx < this->decays.size(); idx++) {
      auto i = this->decays[idx];
      if (i->core != this->top) {
        auto s = i->core->to_string();
        if (data->data_p.find(s) != data->data_p.end()) {
          auto tmp1 = i->core->get_amp(n, data->data_p[s]);
          for (int j = 0; j < n; j++) {
            ret.ptr[j] *= tmp1.ptr[j];
          }
        } else {
          std::cout << s << " not found" << std::endl;
        }
      }
      auto tmp2 = i->get_amp(n, data->data[idx]);
      for (int j = 0; j < n; j++) {
        ret.ptr[j] *= tmp2.ptr[j];
      }
    }
    return ret;
  };

  std::string to_string() {
    std::ostringstream ostr;
    ostr << "[";
    int idx = 0;
    for (auto i : this->decays) {
      if (idx == 0) {
        ostr << i->to_string();
      } else {

        ostr << "," << i->to_string();
      }
      idx++;
    }
    ostr << "]";
    return ostr.str();
  }
};

class DecayChain : public BaseDecayChain {
public:
  DecayChain(std::vector<BaseDecay *> decs) : BaseDecayChain(decs){};
};

class BaseDecayGroup {
public:
  BaseDecayGroup(std::vector<BaseDecayChain *> decs) : decs(decs){};
  std::vector<BaseDecayChain *> decs;
  virtual void init_params(VarManager *vm) {
    for (auto i : this->decs) {
      i->init_params(vm);
    }
  }
  virtual Tensor<std::complex<double>>
  get_amp(size_t n, std::map<std::string, ChainData *> data) {
    Tensor<std::complex<double>> ret(n);
    for (auto i : this->decs) {
      auto s = i->to_string();
      auto it = data.find(s);
      if (it != data.end()) {
        auto tmp = i->get_amp(n, data.at(s));
        for (int j = 0; j < n; j++) {
          ret.ptr[j] += tmp.ptr[j];
        }
      } else {
        std::cout << s << "not found" << std::endl;
      }
    }
    return ret;
  };
  virtual Tensor<double> get_amp2s(size_t n,
                                   std::map<std::string, ChainData *> data) {
    auto ret = Tensor<double>(n);
    auto amp = this->get_amp(n, data);
    for (int i = 0; i < n; i++) {
      ret.ptr[i] = norm(amp.ptr[i]);
    }
    return ret;
  };
};

class DecayGroup : public BaseDecayGroup {
public:
  DecayGroup(std::vector<BaseDecayChain *> decs) : BaseDecayGroup(decs){};
};

#endif
