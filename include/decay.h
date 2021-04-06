#include "data.h"
#include "variable.h"
#include <complex>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#ifndef M_PWADECAY_H_
#define M_PWADECAY_H_ 1

#include "dfunction.h"
class BaseParticle {
public:
  int J, P;
  BaseParticle(std::string name, int J, int P) : name(name), J(J), P(P){};
  BaseParticle(std::string name) : name(name), J(0), P(-1){};
  std::string name;
  virtual void init_params(VarManager *vm) {}
  virtual Tensor<std::complex<double>> get_amp(size_t n, Tensor<double> *data);
  std::string to_string() { return this->name; }
};

const std::function<double()> zeros_f = [] { return 0.0; };

class Particle : public BaseParticle {
public:
  Particle(std::string name, int J = 0, int P = -1)
      : BaseParticle(name, J, P){};
  std::function<double()> mass = zeros_f;
  std::function<double()> width = zeros_f;
  virtual void init_params(VarManager *vm) {
    this->mass = vm->add_var(this->to_string() + "_mass");
    this->width = vm->add_var(this->to_string() + "_width");
  }

  virtual Tensor<std::complex<double>> get_amp(size_t n,
                                               Tensor<double> *data) override;
};

class BaseDecay {
public:
  BaseDecay(BaseParticle *core, std::vector<BaseParticle *> outs)
      : core(core), outs(outs){};
  BaseParticle *core;
  std::vector<BaseParticle *> outs;
  virtual void init_params(VarManager *vm) {}
  virtual SharedTensor get_amp(size_t n, DecayData *data) {
    auto ret = SharedTensor(new ComplexTensor(n));
    for (int i = 0; i < n; i++) {
      ret->ptr[i] = data->angle->alpha->ptr[i];
    }
    return ret;
  };
  std::string to_string() {
    std::ostringstream ostr;
    ostr << this->core->to_string() << "->";
    int idx = 0;
    for (auto i : this->outs) {
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
  bool p_break = false;
  std::vector<std::pair<std::function<double()>, std::function<double()>>> gls;
  Decay(BaseParticle *core, std::vector<BaseParticle *> outs,
        bool p_break = false)
      : BaseDecay(core, outs), p_break(p_break), gls(){};
  virtual void init_params(VarManager *vm) override {
    auto ls = this->get_ls_list();
    this->gls.clear();
    for (int i = 0; i < ls.size(); i++) {
      auto a =
          vm->add_var(this->to_string() + "_g_ls_" + std::to_string(i) + "r");
      auto b =
          vm->add_var(this->to_string() + "_g_ls_" + std::to_string(i) + "i");
      gls.push_back({a, b});
    }
  }

  virtual SharedTensor get_amp(size_t n, DecayData *data) override;
  Tensor<std::complex<double>> get_d_matrix(size_t n, EularAngle *data);

  std::vector<std::pair<int, int>> get_ls_list();

  Tensor<std::complex<double>> get_ls_matrix();

  Tensor<std::complex<double>> get_helicity_amp(size_t n, Tensor<double> *m);
};

class BaseDecayChain {
public:
  BaseDecayChain(std::vector<BaseDecay *> decs) : decays(decs), outs() {
    this->top = decs[0]->core;
    for (auto i : decs) {
      for (auto j : i->outs) {
        bool found = false;
        for (auto k : decs) {
          if (k->core == j)
            found = true;
        }
        if (!found)
          this->outs.push_back(j);
      }
    }
  };
  BaseParticle *top;
  std::vector<BaseParticle *> outs;
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

  virtual Tensor<std::complex<double>> get_amp(size_t n, ChainData *data);
  virtual Tensor<std::complex<double>> get_amp_particle(size_t n,
                                                        ChainData *data);
  virtual Tensor<std::complex<double>> get_amp_decay(size_t n, ChainData *data);
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
  get_amp(size_t n, std::map<std::string, ChainData *> data);

  virtual Tensor<double> get_amp2s(size_t n,
                                   std::map<std::string, ChainData *> data);
};

class DecayGroup : public BaseDecayGroup {
public:
  DecayGroup(std::vector<BaseDecayChain *> decs) : BaseDecayGroup(decs){};
};

#endif
