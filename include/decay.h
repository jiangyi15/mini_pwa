#include "data.h"
#include "variable.h"
#include <complex>
#include <functional>
#include <iostream>
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
  BaseDecay(BaseParticle *core, std::vector<BaseParticle *> outs)
      : core(core), outs(outs){};
  BaseParticle *core;
  std::vector<BaseParticle *> outs;
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
  Decay(BaseParticle *core, std::vector<BaseParticle *> outs)
      : BaseDecay(core, outs){};
  virtual Tensor<std::complex<double>> get_amp(size_t n,
                                               DecayData *data) override {
    auto ret = this->get_d_matrix(n, data->angle);
    auto h = this->get_helicity_amp(n, data->data_p);
    int ja, jb, jc;
    ja = this->core->J;
    jb = this->outs[0]->J;
    jc = this->outs[1]->J;
    for (int i = 0; i < n; i++) {
      for (int ma = -ja; ma <= ja; ma++) {
        for (int mb = -jb; mb <= jb; mb++) {
          for (int mc = -jc; mc <= jc; mc++) {
            if (mb - mc <= ja) {
              ret[i][ma + ja][mb + jb].ptr[mc + jc] *=
                  h[i][mb + jb].ptr[mc + jc];
            }
          }
        }
      }
    }
    return ret;
  };
  Tensor<std::complex<double>> get_d_matrix(size_t n, EularAngle *data) {
    int ja, jb, jc;
    ja = this->core->J;
    jb = this->outs[0]->J;
    jc = this->outs[1]->J;
    auto ret =
        Tensor<std::complex<double>>({n, 2 * ja + 1, 2 * jb + 1, 2 * jc + 1});
    auto w = small_d_weight(2 * ja);
    auto d = small_d_function(2 * ja, &w, data->beta);
    for (int i = 0; i < n; i++) {
      double alpha = data->alpha->ptr[i];
      double gamma = data->gamma->ptr[i];
      for (int ma = -ja; ma <= ja; ma++) {
        for (int mb = -jb; mb <= jb; mb++) {

          for (int mc = -jc; mc <= jc; mc++) {
            if (mb - mc <= ja) {
              ret[i][ma + ja][mb + jb].ptr[mc + jc] =
                  d[i][ma + ja].ptr[mb - mc + ja] *
                  exp(std::complex(0., ma * alpha)) *
                  exp(std::complex(0., (mb - mc) * gamma));
            }
          }
        }
      }
    }
    return ret;
  }

  std::vector<std::pair<int, int>> get_ls_list();

  Tensor<std::complex<double>> get_ls_matrix();

  Tensor<std::complex<double>> get_helicity_amp(size_t n, Tensor<double> *m) {
    int ja, jb, jc;
    ja = this->core->J;
    jb = this->outs[0]->J;
    jc = this->outs[1]->J;
    auto ls_matrix = this->get_ls_matrix();
    auto ret = Tensor<std::complex<double>>({n, 2 * jb + 1, 2 * jc + 1});
    for (int i = 0; i < n; i++) {
      for (int mb = -jb; mb <= jb; mb++) {
        for (int mc = -jc; mc <= jc; mc++) {
          ret[i][mb + jb].ptr[mc + jc] = 0;
          for (int l = 0; l < this->get_ls_list().size(); l++) {
            ret[i][mb + jb].ptr[mc + jc] += ls_matrix[mb + jb][mc + jc].ptr[l];
          }
        }
      }
    }
    return ret;
  }
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
      ret.ptr[i] = norm(amp.ptr[i]) + 1e-10;
    }
    return ret;
  };
};

class DecayGroup : public BaseDecayGroup {
public:
  DecayGroup(std::vector<BaseDecayChain *> decs) : BaseDecayGroup(decs){};
};

#endif
