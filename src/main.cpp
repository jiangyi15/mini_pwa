#include "cnpy.h"
#include "decay.h"
#include "fcn.h"

Tensor<double> *load_array(cnpy::npz_t &npz, std::string idx) {
  cnpy::NpyArray arr = npz[idx];
  double *mv1 = arr.data<double>();
  auto beta = new Tensor<double>(arr.shape, arr.data<double>());
  return beta;
}

std::map<std::string, ChainData *> load_data(cnpy::npz_t &npz) {

  auto alpha = load_array(
      npz, "decay/[A->(B, C)+D, (B, C)->B+C]/A->(B, C)+D/(B, C)/ang/alpha");
  auto beta = load_array(
      npz, "decay/[A->(B, C)+D, (B, C)->B+C]/A->(B, C)+D/(B, C)/ang/beta");
  auto gamma = load_array(
      npz, "decay/[A->(B, C)+D, (B, C)->B+C]/A->(B, C)+D/(B, C)/ang/gamma");

  auto alpha2 = load_array(
      npz, "decay/[A->(B, C)+D, (B, C)->B+C]/(B, C)->B+C/B/ang/alpha");
  auto beta2 = load_array(
      npz, "decay/[A->(B, C)+D, (B, C)->B+C]/(B, C)->B+C/B/ang/beta");
  auto gamma2 = load_array(
      npz, "decay/[A->(B, C)+D, (B, C)->B+C]/(B, C)->B+C/B/ang/gamma");

  auto ma = load_array(npz, "particle/A/m");
  auto mr = load_array(npz, "particle/(B, C)/m");
  auto pa =
      load_array(npz, "decay/[A->(B, C)+D, (B, C)->B+C]/A->(B, C)+D/|q|2");
  auto pr =
      load_array(npz, "decay/[A->(B, C)+D, (B, C)->B+C]/(B, C)->B+C/|q|2");

  auto a1 = new EularAngle(alpha, beta, gamma);
  auto a2 = new EularAngle(alpha2, beta2, gamma2);
  auto d1 = new DecayData(a1, pa);
  auto d2 = new DecayData(a2, pr);

  auto dc1 = new ChainData({d1, d2}, {{"a", ma}, {"s", mr}});

  auto alpha_2 = load_array(
      npz, "decay/[A->(B, D)+C, (B, D)->B+D]/A->(B, D)+C/(B, D)/ang/alpha");
  auto beta_2 = load_array(
      npz, "decay/[A->(B, D)+C, (B, D)->B+D]/A->(B, D)+C/(B, D)/ang/beta");
  auto gamma_2 = load_array(
      npz, "decay/[A->(B, D)+C, (B, D)->B+D]/A->(B, D)+C/(B, D)/ang/gamma");

  auto alpha2_2 = load_array(
      npz, "decay/[A->(B, D)+C, (B, D)->B+D]/(B, D)->B+D/B/ang/alpha");
  auto beta2_2 = load_array(
      npz, "decay/[A->(B, D)+C, (B, D)->B+D]/(B, D)->B+D/B/ang/beta");
  auto gamma2_2 = load_array(
      npz, "decay/[A->(B, D)+C, (B, D)->B+D]/(B, D)->B+D/B/ang/gamma");

  auto ma_2 = load_array(npz, "particle/A/m");
  auto mr_2 = load_array(npz, "particle/(B, D)/m");
  auto pa_2 =
      load_array(npz, "decay/[A->(B, D)+C, (B, D)->B+D]/A->(B, D)+C/|q|2");
  auto pr_2 =
      load_array(npz, "decay/[A->(B, D)+C, (B, D)->B+D]/(B, D)->B+D/|q|2");

  auto a1_2 = new EularAngle(alpha_2, beta_2, gamma_2);
  auto a2_2 = new EularAngle(alpha2_2, beta2_2, gamma2_2);
  auto d1_2 = new DecayData(a1_2, pa_2);
  auto d2_2 = new DecayData(a2_2, pr_2);

  auto dc2 = new ChainData({d1_2, d2_2}, {{"a", ma_2}, {"r", mr_2}});

  auto data = std::map<std::string, ChainData *>(
      {{"[a->s+d,s->b+c]", dc1}, {"[a->r+c,r->b+d]", dc2}});
  return data;
}

int main() {
  cnpy::npz_t npz1 = cnpy::npz_load("../examples/toy.npz");
  auto toy = load_data(npz1);

  cnpy::npz_t npz2 = cnpy::npz_load("../examples/phsp.npz");
  auto phsp = load_data(npz2);

  auto a = new Particle("a");
  auto b = new Particle("b");
  auto c = new Particle("c");
  auto d = new Particle("d");
  auto r = new Particle("r", 1, -1);
  auto r2 = new Particle("s", 1, -1);

  auto dec1 = new Decay(a, {r2, d});
  auto dec2 = new Decay(r2, {b, c});
  auto decs = new DecayChain({dec1, dec2});
  auto dec3 = new Decay(a, {r, c});
  auto dec4 = new Decay(r, {b, d});
  auto decs2 = new DecayChain({dec3, dec4});

  auto dg = new DecayGroup({decs, decs2});

  AmplitudeModel amp(dg);

  amp.fix_params["[a->s+d,s->b+c]total_i"] = 0.0;
  amp.fix_params["[a->s+d,s->b+c]total_r"] = 1.593208560457262;
  amp.fix_params["a->r+c_g_ls_0i"] = 0.0;
  amp.fix_params["a->r+c_g_ls_0r"] = 1.0;
  amp.fix_params["a->s+d_g_ls_0i"] = 0.0;
  amp.fix_params["a->s+d_g_ls_0r"] = 1.0;
  amp.fix_params["r->b+d_g_ls_0i"] = 0.0;
  amp.fix_params["r->b+d_g_ls_0r"] = 1.0;
  amp.fix_params["s->b+c_g_ls_0i"] = 0.0;
  amp.fix_params["s->b+c_g_ls_0r"] = 1.0;
  amp.fix_params["r_mass"] = 1.0;
  amp.fix_params["s_mass"] = 1.0;
  amp.fix_params["r_width"] = 0.1;
  amp.bound["s_width"] = {0.001, 0.2};
  // amp.bound["r_width"] = {0.001, 0.2};
  amp.vm.set_params(amp.fix_params);
  amp.vm.set_params({
      {"r_width", 0.1},
      {"s_width", 0.05},
      {"[a->r+c,r->b+d]total_i", 0.6649200634780428},
      {"[a->r+c,r->b+d]total_r", 0.1193709361286075},
  });

  auto n_data = npz1["particle/C/m"].shape[0];
  auto n_phsp = npz2["particle/C/m"].shape[0];
  std::cout << "n_data: " << n_data << " n_phsp: " << n_phsp << std::endl;

  amp.set_data(n_data, toy);
  amp.set_phsp(n_phsp, phsp);

  amp.vm.save_params("init_params.json");
  auto ai = amp.decay_group->get_amp2s(npz1["particle/C/m"].shape[0], toy);
  for (int i = 0; i < 10; i++) {
    std::cout << ai.ptr[i] << ", ";
  }
  std::cout << std::endl;
  auto min = amp.minimize();
  std::cout << min << std::endl;
  return 0;
}
