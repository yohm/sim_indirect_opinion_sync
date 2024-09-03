#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <queue>
#include <utility>
#include <string>
#include <omp.h>
#include <nlohmann/json.hpp>
#include "PrivRepGame.hpp"


void PrintInitialTimeSeries(PrivateRepGame &prg, const nlohmann::json &params, std::ostream &out = std::cout) {
  size_t num_prints = 20;
  size_t interval = params.at("t_init").get<size_t>() / num_prints;

  size_t N_strategies = prg.Population().size();  // number of different strategies

  for (size_t t = 0; t < num_prints; t++) {
    prg.Update(interval, params.at("q"), params.at("mu_impl"), params.at("mu_assess"), params.at("tau"), false);
    size_t time = (t + 1) * interval;
    out << time << ' ' << prg.SystemWideCooperationLevel();
    if (N_strategies > 1) {
      auto c_levels = prg.NormCooperationLevels();
      for (size_t i = 0; i < N_strategies; i++) {
        for (size_t j = 0; j < N_strategies; j++) {
          out << ' ' << c_levels[i][j];
        }
      }
    }
    out << std::endl;
  }

  prg.ResetCounts();
}

struct OutputMonomorphic {
  double h;
  double hG;
  double h_expected;
  double bc_lower_theory;
  double bc_upper_theory;
};

// monomorphic simulation
OutputMonomorphic RunMonomorphicSimulation(const Norm &norm, const nlohmann::json &params, bool verbose = false) {
  OutputMonomorphic output{};
  PrivateRepGame::population_t pop;
  pop.emplace_back(norm, params["N"].get<size_t>());

  PrivateRepGame prg(pop, params["seed"].get<uint64_t>());
  if (verbose) {
    PrintInitialTimeSeries(prg, params, std::cerr);
  } else {
    prg.Update(params.at("t_init"), params.at("q"), params.at("mu_impl"), params.at("mu_assess"), params.at("tau"),
               false);
  }
  prg.ResetCounts();

  // calculate goodness for homogeneous population
  bool count_goodness = true;
  prg.Update(params["t_measure"].get<size_t>(), params["q"].get<double>(), params["mu_impl"].get<double>(),
             params["mu_assess"].get<double>(), params["tau"].get<double>(), count_goodness);
  // std::cerr << "pc_res_res: " << prg.SystemWideCooperationLevel() << std::endl;

  auto g_g2_vec = prg.NormAverageGoodness();
  assert(g_g2_vec.size() == 1);
  double g = g_g2_vec[0].first;  // == h
  double g2 = g_g2_vec[0].second;
  auto N = static_cast<double>(pop.at(0).second);
  double h = g;
  double hG = (N - 1.0) / (N - 2.0) * g2 / g - 1.0 / (N - 2.0);
  // std::cerr << "h: " << g << std::endl << "hG: " << hG << "\n";
  output.h = g;
  output.hG = hG;

  // expected results
  // h = h hG R_P(G,G) + h(1-hG) (R_P(G,B) + R_P(B,G)) + (1-2h+h hG) R_P(B,B)
  constexpr Reputation G = Reputation::G, B = Reputation::B;
  constexpr Action C = Action::C, D = Action::D;
  const Norm rescaled = norm.RescaleWithError(0.0, params["mu_assess"].get<double>());
  // ActionRule is DISC
  double RP_GG = rescaled.Gprob(G, C);
  double RP_GB = rescaled.Gprob(G, D);
  double RP_BG = rescaled.Gprob(B, C);
  double RP_BB = rescaled.Gprob(B, D);
  double h_expected = h * hG * RP_GG + h * (1.0 - hG) * (RP_GB + RP_BG) + (1.0 - 2.0 * h + h * hG) * RP_BB;
  // std::cerr << "h_expected: " << h_expected << std::endl;
  output.h_expected = h_expected;

  // Delta_G = [R(G,C) - R(G,D)][P(C)-P(D)]
  // Delta_B = [R(B,C) - R(B,D)][P(C)-P(D)]
  double delta_G = rescaled.Gprob(G, C) - rescaled.Gprob(G, D);
  double delta_B = rescaled.Gprob(B, C) - rescaled.Gprob(B, D);
  // alpha_G > 0  <=> [hG delta_G + (1-hG) delta_B]b > c
  double x = hG * delta_G + (1.0 - hG) * delta_B;
  double bc_lower = (x > 0) ? 1.0 / x : std::numeric_limits<double>::infinity();  // if x<0, never stable

  // alpha_B < 0  <=>  [ h(1-hG) delta_G + (1-2h+hhG) delta_B ]b < (1-h)c
  double y = h * (1.0 - hG) * delta_G + (1.0 - 2.0 * h + h * hG) * delta_B;
  double bc_upper = (y > 0) ? (1.0 - h) / y : std::numeric_limits<double>::infinity();  // if y<0, always stable

  // std::cerr << "bc_lower_theory: " << bc_lower << std::endl;
  // std::cerr << "bc_upper_theory: " << bc_upper << std::endl;
  output.bc_lower_theory = bc_lower;
  output.bc_upper_theory = bc_upper;

  return output;
}

struct PolymorphicOutput {
  double bc_lower;
  double bc_upper;
};

// polymorphic simulation with ALLC/ALLD mutants
PolymorphicOutput
RunPolymorphicSimulation(const Norm &norm, const nlohmann::json &params, const Norm &mutant, bool verbose = false) {

  double bc_lower = 1.0;
  double bc_upper = std::numeric_limits<double>::infinity();

  PrivateRepGame::population_t pop;
  pop.emplace_back(norm, params["N"].get<size_t>() - 1);
  pop.emplace_back(mutant, 1);

  PrivateRepGame prg(pop, params["seed"].get<uint64_t>());
  prg.Update(params.at("t_init"), params.at("q"), params.at("mu_impl"), params.at("mu_assess"), params.at("tau"),
             false);
  prg.ResetCounts();

  prg.Update(params["t_measure"].get<size_t>(), params["q"].get<double>(), params["mu_impl"].get<double>(),
             params["mu_assess"].get<double>(), params["tau"].get<double>(), false);
  auto pcs = prg.NormCooperationLevels();
  double pc_res_res = pcs[0][0];
  double pc_res_mut = pcs[0][1];
  double pc_mut_res = pcs[1][0];
  if (verbose) {
    std::cerr << "  pc_res_res: " << pc_res_res << std::endl;
    std::cerr << "  pc_res_mut: " << pc_res_mut << std::endl;
    std::cerr << "  pc_mut_res: " << pc_mut_res << std::endl;
  }
  // resident is stable if b(pc_res_res - pc_res_mut) > c(pc_res_res - pc_mut_res)

  if (pc_res_res - pc_res_mut > 0) {
    double bc_lower_mut = (pc_res_res - pc_mut_res) / (pc_res_res - pc_res_mut);
    if (verbose) std::cerr << "  bc_lower_mut: " << bc_lower_mut << std::endl;
    bc_lower = std::max(bc_lower, bc_lower_mut);
  } else if (pc_res_res - pc_res_mut == 0) {
    if (pc_res_res - pc_mut_res > 0) {
      if (verbose) std::cerr << "  bc_lower_mut: infinity" << std::endl;
      bc_lower = std::numeric_limits<double>::infinity();
    } else {
      if (verbose) std::cerr << "  bc_lower_mut: 1" << std::endl;
    }
  } else {
    double bc_upper_mut = (pc_res_res - pc_mut_res) / (pc_res_res - pc_res_mut);
    if (verbose) std::cerr << "  bc_upper_mut: " << bc_upper_mut << std::endl;
    bc_upper = std::min(bc_upper, bc_upper_mut);
  }

  PolymorphicOutput output;
  output.bc_lower = bc_lower;
  output.bc_upper = bc_upper;

  return output;
}

std::pair<double,double> AverageAndError(const std::vector<double>& vals) {
  double sum = 0.0;
  double n = static_cast<double>(vals.size());
  for (double val : vals) {
    sum += val;
  }
  double avg = sum / n;
  double err = 0.0;
  for (double val : vals) {
    err += (val - avg) * (val - avg);
  }
  if (n <= 1.0) return {avg, std::numeric_limits<double>::infinity()};
  err = std::sqrt(err / (n * (n-1)));
  return {avg, err};
}

int main(int argc, char *argv[]) {

  std::queue<std::string> args;
  nlohmann::json params = nlohmann::json::object();
  // -j param.json : set parameters by json file or json string
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "-j" && i + 1 < argc) {
      std::ifstream fin(argv[++i]);
      // check if file exists
      if (fin) {
        fin >> params;
        fin.close();
      } else {
        std::istringstream iss(argv[i]);
        iss >> params;
      }
    } else {
      args.emplace(argv[i]);
    }
  }

  // set default parameters
  const nlohmann::json default_params = {{"N",         100},
                                         {"t_init",    1e3},
                                         {"t_measure", 1e3},
                                         {"q",         1.0},
                                         {"mu_impl",   0.0},
                                         {"mu_assess", 0.02},
                                         {"tau",       0.0},
                                         {"seed",      123456789ull}};
  for (auto it = default_params.begin(); it != default_params.end(); ++it) {
    if (!params.contains(it.key())) {
      params[it.key()] = it.value();
    }
  }

  // if params has keys other than default_params, print them
  for (auto it = params.begin(); it != params.end(); ++it) {
    if (!default_params.contains(it.key())) {
      std::cerr << "[Error] unknown parameter: " << it.key() << std::endl;
      return 1;
    }
  }

  std::cerr << "Parameters:" << std::endl;
  std::cerr << params.dump(2) << std::endl;

  std::cerr << "OMP_NUM_THREADS: " << omp_get_max_threads() << std::endl;

  auto show_usage = [&argv, default_params] {
    std::cerr << "Usage: " << argv[0] << " [-j param.json] norm1 num_samples" << std::endl;
    std::cerr << "Default parameters:" << std::endl;
    std::cerr << "  " << default_params.dump(2) << std::endl;
  };


  if (args.size() != 2) {
    std::cerr << "[Error] wrong input format" << std::endl;
    show_usage();
    return 1;
  }
  Norm norm = Norm::ParseNormString(args.front(), false);
  size_t num_samples = std::stoul(args.back());

  if (norm.P.CProb(Reputation::G) != 1.0 || norm.P.CProb(Reputation::B) != 0.0) {
    std::cerr << "[Error] action rule must be the discriminator (1,0)" << std::endl;
    return 1;
  }
  if (params.at("mu_impl").get<double>() != 0.0) {
    std::cerr << "[Error] non-zero mu_impl is not implemented yet" << std::endl;
    return 1;
  }

  auto start = std::chrono::high_resolution_clock::now();

  std::vector<OutputMonomorphic> outputs(num_samples);
  #pragma omp parallel for
  for (size_t i = 0; i < num_samples; i++) {
    nlohmann::json params_copy = params;  // copy parameter
    params_copy["seed"] = params["seed"].get<uint64_t>() + i;
    auto output = RunMonomorphicSimulation(norm, params_copy, false);
    outputs[i] = output;
  }

  std::vector<PolymorphicOutput> outputs_alld(num_samples);
  #pragma omp parallel for
  for (size_t i = 0; i < num_samples; i++) {
    nlohmann::json params_copy = params;  // copy parameter
    params_copy["seed"] = params["seed"].get<uint64_t>() + i + 123456789;
    auto out_alld = RunPolymorphicSimulation(norm, params_copy, Norm::AllD(), false);
    outputs_alld[i] = out_alld;
  }

  std::vector<PolymorphicOutput> outputs_allc(num_samples);
  #pragma omp parallel for
  for (size_t i = 0; i < num_samples; i++) {
    nlohmann::json params_copy = params;  // copy parameter
    params_copy["seed"] = params["seed"].get<uint64_t>() + i + 987654321;
    auto out_allc = RunPolymorphicSimulation(norm, params_copy, Norm::AllC(), false);
    outputs_allc[i] = out_allc;
  }

  nlohmann::json results;
  {
    std::vector<double> h, hG, h_expected, bc_lower_theory, bc_upper_theory;
    for (const auto &output : outputs) {
      h.push_back(output.h);
      hG.push_back(output.hG);
      h_expected.push_back(output.h_expected);
      bc_lower_theory.push_back(output.bc_lower_theory);
      bc_upper_theory.push_back(output.bc_upper_theory);
    }
    auto [h_avg, h_err] = AverageAndError(h);
    auto [hG_avg, hG_err] = AverageAndError(hG);
    auto [h_expected_avg, h_expected_err] = AverageAndError(h_expected);
    auto [bc_lower_theory_avg, bc_lower_theory_err] = AverageAndError(bc_lower_theory);
    auto [bc_upper_theory_avg, bc_upper_theory_err] = AverageAndError(bc_upper_theory);
    results = {{"h",               h_avg},
               {"h_err",           h_err},
               {"hG",              hG_avg},
               {"hG_err",          hG_err},
               {"h_expected",      h_expected_avg},
               {"h_expected_err",  h_expected_err},
               {"bc_lower_theory", bc_lower_theory_avg},
               {"bc_lower_theory_err", bc_lower_theory_err},
               {"bc_upper_theory", bc_upper_theory_avg},
               {"bc_upper_theory_err", bc_upper_theory_err}};
  }
  {
    std::vector<double> bc_lower;
    for (const auto &output : outputs_alld) {
      bc_lower.push_back(output.bc_lower);
    }
    auto [bc_lower_avg, bc_lower_err] = AverageAndError(bc_lower);
    results["bc_lower_alld"] = bc_lower_avg;
    results["bc_lower_alld_err"] = bc_lower_err;
  }
  {
    std::vector<double> bc_upper;
    for (const auto &output : outputs_allc) {
      bc_upper.push_back(output.bc_upper);
    }
    auto [bc_upper_avg, bc_upper_err] = AverageAndError(bc_upper);
    results["bc_upper_allc"] = bc_upper_avg;
    results["bc_upper_allc_err"] = bc_upper_err;
  }

  std::cout << results.dump(2) << std::endl;

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cerr << "Elapsed time: " << elapsed.count() << " s\n";

  return 0;
}