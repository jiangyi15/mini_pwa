#include "decay.h"

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

std::vector<std::pair<int, int>> Decay::get_ls_list()

{
  int ja, jb, jc;
  ja = this->core->J;
  jb = this->outs[0]->J;
  jc = this->outs[1]->J;
  int pa, pb, pc;
  pa = this->core->P;
  pb = this->outs[0]->P;
  pc = this->outs[1]->P;
  return ::get_ls_list(ja, pa, jb, pb, jc, pc);
}

Tensor<std::complex<double>> Decay::get_ls_matrix() {
  int ja, jb, jc;
  ja = this->core->J;
  jb = this->outs[0]->J;
  jc = this->outs[1]->J;
  auto ls = this->get_ls_list();
  auto ret = Tensor<std::complex<double>>({2 * jb + 1, 2 * jc + 1, ls.size()});
  return ret;
}
