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
  virtual std::complex<double> get_amp(Tensor<double> *m) { return 1.0; };
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
  std::complex<double> get_amp(Tensor<double> *data) {
    double m = data->ptr[0];
    std::cout << m << std::endl;
    double m0 = this->mass();
    double g0 = this->width();
    return 1. / std::complex(m0 * m0 - m * m, -m0 * g0);
  };
};

class BaseDecay {
public:
  BaseDecay(BaseParticle *core, std::vector<BaseParticle *> decay)
      : core(core), decay(decay){};
  BaseParticle *core;
  std::vector<BaseParticle *> decay;
  virtual void init_params(VarManager *vm) {}
  virtual std::complex<double> get_amp(DecayData *data) {
    return data->angle->alpha->ptr[0];
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
  BaseDecayChain(std::vector<BaseDecay *> decs) : decays(decs){};
  std::vector<BaseDecay *> decays;
  virtual void init_params(VarManager *vm) {
    for (auto i : this->decays) {
      i->init_params(vm);
      i->core->init_params(vm);
    }
  }

  virtual std::complex<double> get_amp(ChainData *data) {
    std::complex<double> a = 1.0;
    int idx = 0;
    for (auto i : this->decays) {
      std::cout << "a" << std::endl;
      a *= i->core->get_amp(data->data_p[i->core->to_string()]) *
           i->get_amp(data->data[idx]);
      std::cout << "a" << std::endl;
      idx += 1;
    }
    return a;
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
  virtual std::complex<double>
  get_amp(std::map<std::string, ChainData *> data) {
    std::complex<double> a = 0.0;
    for (auto i : this->decs) {
      auto s = i->to_string();
      a += i->get_amp(data.at(s));
    }
    return a;
  };
  virtual double get_amp2s(std::map<std::string, ChainData *> data) {
    auto a = abs(this->get_amp(data));
    return a * a;
  };
};

class DecayGroup : public BaseDecayGroup {
public:
  DecayGroup(std::vector<BaseDecayChain *> decs) : BaseDecayGroup(decs){};
};

#endif
