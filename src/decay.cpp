#include "decay.h"
#include "cg_coeff.h"

Tensor<std::complex<double>> BaseParticle::get_amp(size_t n, Tensor<double> *m,
                                                   Tensor<double> *q) {
  Tensor<std::complex<double>> ret(n);
  for (int i = 0; i < n; i++) {
    ret.ptr[i] = 1.;
  }
  return ret;
};

Tensor<std::complex<double>> Particle::get_amp(size_t n, Tensor<double> *m,
                                               Tensor<double> *q) {
  Tensor<std::complex<double>> ret(n);
  double m0 = this->mass();
  double g0 = this->width();
  for (int i = 0; i < n; i++) {
    double mi = m->ptr[i];
    // if (i==0) std::cout << this->name <<" "<< mi <<" " << m0 << " " << g0 <<
    // " ";
    ret.ptr[i] = 1. / std::complex<double>(m0 * m0 - mi * mi, -m0 * g0);
  }
  // std::cout << ret << " "<< *m << " " << m0 << " " << g0 << std::endl;
  // exit(0);
  return ret;
};

std::vector<std::pair<int, int>> get_ls_list(int ja, int pa, int jb, int pb,
                                             int jc, int pc,
                                             bool p_break = false) {
  std::vector<std::pair<int, int>> ret;
  auto dl = 0;
  if (!p_break) {
    if ((pa * pb * pc) == 1)
      dl = 0;
    else
      dl = 1;
  } //  pa = pb * pc * (-1)^l
  auto s_min = abs(jb - jc);
  auto s_max = jb + jc;
  for (int s = s_min; s <= s_max; s++) {
    for (int l = abs(ja - s); l <= ja + s; l++) {
      if (!p_break) {
        if (l % 2 == dl) {
          ret.push_back({l, s});
        }
      } else {
        ret.push_back({l, s});
      }
    }
  }
  return ret;
}

std::vector<std::pair<int, int>> Decay::get_ls_list() {
  int ja, jb, jc;
  ja = this->core->J;
  jb = this->outs[0]->J;
  jc = this->outs[1]->J;
  int pa, pb, pc;
  pa = this->core->P;
  pb = this->outs[0]->P;
  pc = this->outs[1]->P;
  auto ret = ::get_ls_list(ja, pa, jb, pb, jc, pc, this->p_break);
  return ret;
}

Tensor<std::complex<double>> Decay::get_ls_matrix() {
  int ja, jb, jc;
  ja = this->core->J;
  jb = this->outs[0]->J;
  jc = this->outs[1]->J;
  auto ls = this->get_ls_list();
  auto ret = Tensor<std::complex<double>>(
      {2 * (size_t)jb + 1, 2 * (size_t)jc + 1, ls.size()});
  for (int i = 0; i < ls.size(); i++) {
    auto l = ls[i].first;
    auto s = ls[i].second;
    for (int lambda_b = -jb; lambda_b <= jb; lambda_b++) {
      for (int lambda_c = -jc; lambda_c <= jc; lambda_c++) {
        ret[lambda_b + jb][lambda_c + jc].ptr[i] =
            sqrt((2.0 * l + 1) / (2.0 * ja + 1)) *
            cg_coeff2(2 * jb, 2 * jc, 2 * lambda_b, -2 * lambda_c, 2 * s,
                      2 * lambda_b - 2 * lambda_c) *
            cg_coeff2(2 * l, 2 * s, 0, 2 * lambda_b - 2 * lambda_c, 2 * ja,
                      2 * lambda_b - 2 * lambda_c);
      }
    }
  }
  return ret;
}

Tensor<std::complex<double>> Decay::get_d_matrix(size_t n, EularAngle *data) {
  int ja, jb, jc;
  ja = this->core->J;
  jb = this->outs[0]->J;
  jc = this->outs[1]->J;
  auto ret = Tensor<std::complex<double>>(
      {n, 2 * (size_t)ja + 1, 2 * (size_t)jb + 1, 2 * (size_t)jc + 1});
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
                exp(std::complex<double>(0., ma * alpha)) *
                exp(std::complex<double>(0., (mb - mc) * gamma));
          }
        }
      }
    }
  }
  return ret;
}

SharedTensor Decay::get_amp(size_t n, DecayData *data, Tensor<double> *m) {
  auto d = this->get_d_matrix(n, data->angle);
  auto h = this->get_helicity_amp(n, m, data->data_p);
  int ja, jb, jc;
  ja = this->core->J;
  jb = this->outs[0]->J;
  jc = this->outs[1]->J;
  auto ret = SharedTensor(new Tensor<std::complex<double>>(
      {n, 2 * (size_t)ja + 1, 2 * (size_t)jb + 1, 2 * (size_t)jc + 1}));
  for (int i = 0; i < n; i++) {
    for (int ma = -ja; ma <= ja; ma++) {
      for (int mb = -jb; mb <= jb; mb++) {
        for (int mc = -jc; mc <= jc; mc++) {
          if (mb - mc <= ja) {
            (*ret)[i][ma + ja][mb + jb].ptr[mc + jc] =
                d[i][ma + ja][mb + jb].ptr[mc + jc] *
                h[i][mb + jb].ptr[mc + jc];
          }
        }
      }
    }
  }
  return ret;
};

Tensor<std::complex<double>>
Decay::get_helicity_amp(size_t n, Tensor<double> *m, Tensor<double> *q2) {
  int ja, jb, jc;
  ja = this->core->J;
  jb = this->outs[0]->J;
  jc = this->outs[1]->J;
  auto ls_matrix = this->get_ls_matrix();
  auto ret =
      Tensor<std::complex<double>>({n, 2 * (size_t)jb + 1, 2 * (size_t)jc + 1});
  auto ls_list = this->get_ls_list();
  for (int i = 0; i < n; i++) {
    for (int mb = -jb; mb <= jb; mb++) {
      for (int mc = -jc; mc <= jc; mc++) {
        ret[i][mb + jb].ptr[mc + jc] = 0.;
        for (int l = 0; l < ls_list.size(); l++) {
          ret[i][mb + jb].ptr[mc + jc] +=
              std::complex<double>(this->gls[l].first(),
                                   this->gls[l].second()) *
              ls_matrix[mb + jb][mc + jc].ptr[l] *
              pow(sqrt(q2->ptr[i]), ls_list[i].first);
        }
      }
    }
  }
  return ret;
}

template <typename K, typename V> void print_map(std::map<K, V> &m) {
  for (auto i : m) {
    std::cout << i.first << ":" << i.second << " ";
  }
  std::cout << std::endl;
}

void small_transpose(std::shared_ptr<ComplexTensor> a, size_t idx1,
                     size_t idx2) {
  /*
   * a[...,idx1,...,idx2,...] <=> a[..., idx2, ...., idx1, ...]
   *    n1  m1   n2  m2  n3          n1   m2    n2    m1   n3
   */
  if (idx1 == idx2)
    return;
  int n1 = 1, n2 = 1, n3 = 1;
  for (int i = 0; i < a->shape.size(); i++) {
    if (i < std::min(idx1, idx2)) {
      n1 *= a->shape[i];
    } else if (i > std::max(idx1, idx2)) {
      n3 *= a->shape[i];
    } else if (i != idx1 && i != idx2) {
      n2 *= a->shape[i];
    }
  }
  auto m1 = a->shape[idx1], m2 = a->shape[idx2];

  for (int i1 = 0; i1 < n1; i1++) {
    for (int i2 = 0; i2 < n2; i2++) {
      for (int i3 = 0; i3 < n3; i3++) {
        for (int j1 = 0; j1 < m1; j1++) {
          for (int j2; j2 < m2; j2++) {
            auto p1 = i3 + n3 * (j2 + m2 * (i2 + n2 * (j1 + m1 * i1)));
            auto p2 = i3 + n3 * (j2 + m1 * (i2 + n2 * (j1 + m2 * i1)));
            auto tmp = a->ptr[p1];
            a->ptr[p1] = a->ptr[p2];
            a->ptr[p2] = tmp;
          }
        }
      }
    }
  }
}

std::shared_ptr<NamedTensor> combine_amp(size_t n, NamedTensor &a,
                                         NamedTensor &b) {
  /**
   * Decay Chain Sum
   * a->r+d, r->b+c
   *
   * A = A(lambda_a, extra1, lambda_r, extra2)
   * B = B(lambda_r, extra3)
   *
   * ret -> sum_r A_a A_r (lambda_a, extra1, extra2)
   *
   * loop: B, lambda_r, extra3
   *       A, lambda_a + extra1, lambda_r, extra2
   * => sum(  r*n_e3 + e3 , i * (nr*n_e2) +  r*n_e2 + n_e2)
   *
   */

  std::map<std::string, size_t, std::less<std::string>> all_shape;

  // print_map(a.index);
  // print_map(b.index);
  for (auto i : a.index) {
    all_shape[i.first] = a.data->shape[i.second];
  }
  std::string sum_index;
  for (auto i : b.index) {
    if (all_shape.find(i.first) != all_shape.end())
      all_shape[i.first] =
          std::max(b.data->shape[i.second], all_shape.at(i.first));
    else
      all_shape[i.first] = b.data->shape[i.second];
    if (i.second == 1)
      sum_index = i.first;
  }

  // print_map(all_shape);

  std::vector<size_t> return_shape({n});
  std::map<std::string, size_t, std::less<std::string>> return_index;
  int idx = 1;
  for (auto i : all_shape) {
    if (i.first != sum_index) {
      return_shape.push_back(i.second);
      return_index[i.first] = i.second;
      idx++;
    }
  }

  int ne1 = 1, ne2 = 1, nr, ne3 = 1;
  nr = all_shape[sum_index];
  for (auto i : b.index) {
    if (i.first != sum_index) {
      ne3 *= b.data->shape[i.second];
    }
  }
  auto a_index = b.index[sum_index];

  for (auto i : a.index) {
    if (i.first != sum_index) {
      if (i.second < a_index)
        ne1 *= a.data->shape[i.second];
      else
        ne2 *= a.data->shape[i.second];
    }
  }

  auto ret = SharedTensor(new Tensor<std::complex<double>>(return_shape));

  for (int idx = 0; idx < n; idx++) {
    for (int e1 = 0; e1 < ne1; e1++) {
      for (int e2 = 0; e2 < ne2; e2++) {
        for (int e3 = 0; e3 < ne3; e3++) {
          int idx_ret = e1 * ne3 * ne2 + e2 * ne3 + e3;
          (*ret)[idx].ptr[idx_ret] = 0.;
          for (int r = 0; r < nr; r++) {
            int idx_a = e1 * nr * ne2 + r * ne2 + e2;
            int idx_b = r * ne3 + e3;
            (*ret)[idx].ptr[idx_ret] +=
                (*a.data)[idx].ptr[idx_a] * (*b.data)[idx].ptr[idx_b];
          }
        }
      }
    }
  }

  // transpose to real_shape;
  std::map<size_t, size_t, std::less<size_t>> trans_index;
  std::map<size_t, size_t, std::less<size_t>> shape_index;
  size_t ptr = 1;
  for (auto i : a.index) {
    if (i.first != sum_index) {
      if (i.second < a.index[sum_index]) {
        auto tmp = 1;
        for (auto j : return_index) {
          if (j.first == i.first) {
            trans_index[ptr] = tmp;
            shape_index[ptr++] = all_shape[i.first];
          }
          tmp++;
        }
      }
    }
  }
  for (auto i : a.index) {
    if (i.first != sum_index) {
      if (i.second > a.index[sum_index]) {
        auto tmp = 1;
        for (auto j : return_index) {
          if (j.first == i.first) {
            trans_index[ptr] = tmp;
            shape_index[ptr++] = all_shape[i.first];
          }
          tmp++;
        }
      }
    }
  }

  size_t ptr_b = ptr;
  for (auto i : b.index) {
    if (i.first != sum_index) {
      auto tmp = 1;
      for (auto j : return_index) {
        if (j.first == i.first) {
          trans_index[ptr] = tmp;
          shape_index[ptr++] = all_shape[i.first];
        }
        tmp++;
      }
    }
  }

  std::map<size_t, size_t, std::less<size_t>> trans_index_need;
  for (auto i : trans_index) {
    if (i.first != i.second) {
      if (shape_index[i.first] != 1) {
        trans_index_need[i.first] = i.second;
      }
    }
  }

  for (auto i : trans_index_need) {
    small_transpose(ret, i.first, i.second);
  }

  return std::shared_ptr<NamedTensor>(new NamedTensor(return_index, ret));
}

Tensor<std::complex<double>> BaseDecayChain::get_amp_particle(size_t n,
                                                              ChainData *data) {
  Tensor<std::complex<double>> ret(n);
  for (int i = 0; i < n; i++) {
    ret.ptr[i] = std::complex<double>(this->total_r(), this->total_i());
  }

  for (int idx = 0; idx < this->decays.size(); idx++) {
    auto i = this->decays[idx];
    auto s = i->core->to_string();
    if (i->core != this->top) {
      if (data->data_p.find(s) != data->data_p.end()) {
        auto tmp1 =
            i->core->get_amp(n, data->data_p[s], data->data[idx]->data_p);
        for (int j = 0; j < n; j++) {
          ret.ptr[j] *= tmp1.ptr[j];
        }
      } else {
        std::cout << s << " not found" << std::endl;
      }
    }
  }
  return ret;
}

Tensor<std::complex<double>> BaseDecayChain::get_amp_decay(size_t n,
                                                           ChainData *data) {

  std::map<std::string, std::shared_ptr<NamedTensor>> decay_amp;
  std::map<std::string, std::map<std::string, int>> index;
  for (int idx = 0; idx < this->decays.size(); idx++) {
    auto i = this->decays[idx];
    auto s = i->core->to_string();
    auto tmp = i->get_amp(n, data->data[idx], data->data_p[s]);
    auto name = std::map<std::string, size_t>{{i->core->to_string(), 1}};
    for (int pi = 0; pi < i->outs.size(); pi++) {
      name[i->outs[pi]->to_string()] = pi + 2;
    }
    decay_amp[s] = std::shared_ptr<NamedTensor>(new NamedTensor(name, tmp));
  }
  // need to keep the order(leaf first)
  for (int idx = this->decays.size() - 1; idx >= 0; idx--) {
    auto i = this->decays[idx];
    auto tmp2 = decay_amp.at(i->core->to_string());
    for (int j = 0; j < i->outs.size(); j++) {
      auto pi = i->outs[j];
      auto si = pi->to_string();
      auto it = decay_amp.find(si);
      if (it != decay_amp.end()) {
        auto amp_i = decay_amp.at(si);
        auto tmp = combine_amp(n, *tmp2, *amp_i);
        decay_amp[i->core->to_string()] = tmp;
      }
    }
  }
  auto ret2 = decay_amp[this->top->to_string()];
  Tensor<std::complex<double>> ret3(n);
  for (int i = 0; i < n; i++) {
    ret3.ptr[i] = ret2->data->ptr[i];
  }
  return ret3;
};

Tensor<std::complex<double>> BaseDecayChain::get_amp(size_t n,
                                                     ChainData *data) {
  auto amp = this->get_amp_particle(n, data);
  auto amp2 = this->get_amp_decay(n, data);

  std::vector<size_t> shape({n});
  size_t total_n = 1;
  for (auto i : this->outs) {
    auto tmp = 2 * (size_t)(i->J) + 1;
    total_n *= tmp;
    shape.push_back(tmp);
  }

  Tensor<std::complex<double>> ret3(shape);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < total_n; j++) {
      ret3.ptr[i * total_n + j] = amp.ptr[i] * amp2.ptr[i * total_n + j];
    }
  }
  return ret3;
};

Tensor<std::complex<double>>
BaseDecayGroup::get_amp(size_t n, std::map<std::string, ChainData *> data) {
  std::vector<size_t> shape({n});
  for (auto i : this->decs[0]->outs)
    shape.push_back(2 * (size_t)(i->J) + 1);
  Tensor<std::complex<double>> ret(shape);
  for (auto i : this->decs) {
    auto s = i->to_string();
    auto it = data.find(s);
    if (it != data.end()) {
      auto tmp = i->get_amp(n, data.at(s));
      for (int j = 0; j < ret.total_shape(); j++) {
        ret.ptr[j] += tmp.ptr[j];
      }
    } else {
      std::cout << s << "not found" << std::endl;
    }
  }
  return ret;
};

Tensor<double>
BaseDecayGroup::get_amp2s(size_t n, std::map<std::string, ChainData *> data) {
  auto ret = Tensor<double>(n);
  auto amp = this->get_amp(n, data);
  size_t total_n = 1;
  for (auto i : this->decs[0]->outs)
    total_n *= 2 * (size_t)(i->J) + 1;
  for (int i = 0; i < n; i++) {
    ret.ptr[i] = 0.;
    for (int j = 0; j < total_n; j++) {
      ret.ptr[i] += norm(amp.ptr[i * total_n + j]);
    }
    if (std::isnan(ret.ptr[i])) {
      ret.ptr[i] = 1e-10;
    } else {
      ret.ptr[i] = std::max(std::min(ret.ptr[i], 1e10), 1e-10);
    }
  }
  std::cout << ret << std::endl;
  return ret;
};
