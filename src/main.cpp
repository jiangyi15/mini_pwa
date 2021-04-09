#include "cnpy.h"
#include "decay.h"
#include "fcn.h"

const std::vector<std::string> all_keys = {
    "decay/[A->(B, D)+C, (B, D)->B+D]/(B, D)->B+D/D/ang/beta",
    "decay/[A->(B, D)+C, (B, D)->B+D]/(B, D)->B+D/D/ang/alpha",
    "decay/[A->(B, D)+C, (B, D)->B+D]/(B, D)->B+D/D/ang/gamma",
    "decay/[A->(B, D)+C, (B, D)->B+D]/(B, D)->B+D/D/aligned_angle/beta",
    "decay/[A->(B, D)+C, (B, D)->B+D]/(B, D)->B+D/D/aligned_angle/alpha",
    "decay/[A->(B, D)+C, (B, D)->B+D]/(B, D)->B+D/D/aligned_angle/gamma",
    "decay/[A->(B, D)+C, (B, D)->B+D]/(B, D)->B+D/B/ang/beta",
    "decay/[A->(B, D)+C, (B, D)->B+D]/(B, D)->B+D/B/ang/alpha",
    "decay/[A->(B, D)+C, (B, D)->B+D]/(B, D)->B+D/B/ang/gamma",
    "decay/[A->(B, D)+C, (B, D)->B+D]/(B, D)->B+D/B/aligned_angle/beta",
    "decay/[A->(B, D)+C, (B, D)->B+D]/(B, D)->B+D/B/aligned_angle/alpha",
    "decay/[A->(B, D)+C, (B, D)->B+D]/(B, D)->B+D/B/aligned_angle/gamma",
    "decay/[A->(B, D)+C, (B, D)->B+D]/A->(B, D)+C/(B, D)/ang/beta",
    "decay/[A->(B, D)+C, (B, D)->B+D]/A->(B, D)+C/(B, D)/ang/alpha",
    "decay/[A->(B, D)+C, (B, D)->B+D]/A->(B, D)+C/(B, D)/ang/gamma",
    "decay/[A->(B, D)+C, (B, D)->B+D]/A->(B, D)+C/C/ang/beta",
    "decay/[A->(B, D)+C, (B, D)->B+D]/A->(B, D)+C/C/ang/alpha",
    "decay/[A->(B, D)+C, (B, D)->B+D]/A->(B, D)+C/C/ang/gamma",
    "decay/[A->(B, C)+D, (B, C)->B+C]/A->(B, C)+D/(B, C)/ang/beta",
    "decay/[A->(B, C)+D, (B, C)->B+C]/A->(B, C)+D/(B, C)/ang/alpha",
    "decay/[A->(B, C)+D, (B, C)->B+C]/A->(B, C)+D/(B, C)/ang/gamma",
    "decay/[A->(B, C)+D, (B, C)->B+C]/A->(B, C)+D/D/ang/beta",
    "decay/[A->(B, C)+D, (B, C)->B+C]/A->(B, C)+D/D/ang/alpha",
    "decay/[A->(B, C)+D, (B, C)->B+C]/A->(B, C)+D/D/ang/gamma",
    "decay/[A->(B, C)+D, (B, C)->B+C]/(B, C)->B+C/C/ang/beta",
    "decay/[A->(B, C)+D, (B, C)->B+C]/(B, C)->B+C/C/ang/alpha",
    "decay/[A->(B, C)+D, (B, C)->B+C]/(B, C)->B+C/C/ang/gamma",
    "decay/[A->(B, C)+D, (B, C)->B+C]/(B, C)->B+C/C/aligned_angle/beta",
    "decay/[A->(B, C)+D, (B, C)->B+C]/(B, C)->B+C/C/aligned_angle/alpha",
    "decay/[A->(B, C)+D, (B, C)->B+C]/(B, C)->B+C/C/aligned_angle/gamma",
    "decay/[A->(B, C)+D, (B, C)->B+C]/(B, C)->B+C/B/ang/beta",
    "decay/[A->(B, C)+D, (B, C)->B+C]/(B, C)->B+C/B/ang/alpha",
    "decay/[A->(B, C)+D, (B, C)->B+C]/(B, C)->B+C/B/ang/gamma",
    "decay/[A->(B, C)+D, (B, C)->B+C]/(B, C)->B+C/B/aligned_angle/beta",
    "decay/[A->(B, C)+D, (B, C)->B+C]/(B, C)->B+C/B/aligned_angle/alpha",
    "decay/[A->(B, C)+D, (B, C)->B+C]/(B, C)->B+C/B/aligned_angle/gamma",
    "decay/[A->(C, D)+B, (C, D)->C+D]/(C, D)->C+D/C/ang/beta",
    "decay/[A->(C, D)+B, (C, D)->C+D]/(C, D)->C+D/C/ang/alpha",
    "decay/[A->(C, D)+B, (C, D)->C+D]/(C, D)->C+D/C/ang/gamma",
    "decay/[A->(C, D)+B, (C, D)->C+D]/(C, D)->C+D/C/aligned_angle/beta",
    "decay/[A->(C, D)+B, (C, D)->C+D]/(C, D)->C+D/C/aligned_angle/alpha",
    "decay/[A->(C, D)+B, (C, D)->C+D]/(C, D)->C+D/C/aligned_angle/gamma",
    "decay/[A->(C, D)+B, (C, D)->C+D]/(C, D)->C+D/D/ang/beta",
    "decay/[A->(C, D)+B, (C, D)->C+D]/(C, D)->C+D/D/ang/alpha",
    "decay/[A->(C, D)+B, (C, D)->C+D]/(C, D)->C+D/D/ang/gamma",
    "decay/[A->(C, D)+B, (C, D)->C+D]/(C, D)->C+D/D/aligned_angle/beta",
    "decay/[A->(C, D)+B, (C, D)->C+D]/(C, D)->C+D/D/aligned_angle/alpha",
    "decay/[A->(C, D)+B, (C, D)->C+D]/(C, D)->C+D/D/aligned_angle/gamma",
    "decay/[A->(C, D)+B, (C, D)->C+D]/A->(C, D)+B/(C, D)/ang/beta",
    "decay/[A->(C, D)+B, (C, D)->C+D]/A->(C, D)+B/(C, D)/ang/alpha",
    "decay/[A->(C, D)+B, (C, D)->C+D]/A->(C, D)+B/(C, D)/ang/gamma",
    "decay/[A->(C, D)+B, (C, D)->C+D]/A->(C, D)+B/B/ang/beta",
    "decay/[A->(C, D)+B, (C, D)->C+D]/A->(C, D)+B/B/ang/alpha",
    "decay/[A->(C, D)+B, (C, D)->C+D]/A->(C, D)+B/B/ang/gamma",
    "particle/D/p",
    "particle/D/m",
    "particle/(C, D)/p",
    "particle/(C, D)/m",
    "particle/(B, C)/p",
    "particle/(B, C)/m",
    "particle/C/p",
    "particle/C/m",
    "particle/A/p",
    "particle/A/m",
    "particle/(B, D)/p",
    "particle/(B, D)/m",
    "particle/B/p",
    "particle/B/m"};

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

  auto a1 = new EularAngle(alpha, beta, gamma);
  auto a2 = new EularAngle(alpha2, beta2, gamma2);
  auto d1 = new DecayData(a1, alpha);
  auto d2 = new DecayData(a2, alpha2);

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

  auto a1_2 = new EularAngle(alpha_2, beta_2, gamma_2);
  auto a2_2 = new EularAngle(alpha2_2, beta2_2, gamma2_2);
  auto d1_2 = new DecayData(a1_2, alpha_2);
  auto d2_2 = new DecayData(a2_2, alpha2_2);

  auto dc2 = new ChainData({d1_2, d2_2}, {{"a", ma_2}, {"r", mr_2}});

  auto data = std::map<std::string, ChainData *>(
      {{"[a->s+d,s->b+c]", dc1}, {"[a->r+c,r->b+d]", dc2}});
  return data;
}

int main() {
  cnpy::npz_t npz1 = cnpy::npz_load("../toy.npz");
  auto toy = load_data(npz1);

  cnpy::npz_t npz2 = cnpy::npz_load("../phsp.npz");
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

  amp.vm.set_params({{"r_width", 0.1}, {"s_width", 0.05}});
  amp.fix_params["[a->r+c,r->b+d]total_i"] = 1.;
  amp.fix_params["[a->r+c,r->b+d]total_r"] = 1.;
  amp.fix_params["a->r+c_g_ls_0i"] = 1.0;
  amp.fix_params["a->r+c_g_ls_0r"] = 1.0;
  amp.fix_params["a->s+d_g_ls_0i"] = 1.0;
  amp.fix_params["a->s+d_g_ls_0r"] = 1.0;
  amp.fix_params["r->b+d_g_ls_0i"] = 1.0;
  amp.fix_params["r->b+d_g_ls_0r"] = 1.0;
  amp.fix_params["s->b+c_g_ls_0i"] = 1.0;
  amp.fix_params["s->b+c_g_ls_0r"] = 1.0;
  amp.fix_params["r_mass"] = 1.0;
  amp.fix_params["s_mass"] = 1.0;
  amp.bound["s_width"] = {0.001, 0.2};

  amp.set_data(npz1["particle/C/m"].shape[0], toy);
  amp.set_phsp(npz1["particle/C/m"].shape[0], phsp);

  amp.vm.save_params("init_params.json");
  auto min = amp.minimize();
  std::cout << min << std::endl;
  return 0;
}
