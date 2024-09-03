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


struct CooperationProbs {
  double pc_res_res;
  double pc_res_mut;
  double pc_mut_res;
  double pc_mut_mut;
};

CooperationProbs RunPolymorphicSimulation(const Norm &resident, const Norm &mutant, size_t num_mutants, const nlohmann::json &params, bool verbose = false) {

  PrivateRepGame::population_t pop;
  pop.emplace_back(resident, params["N"].get<size_t>() - num_mutants);
  pop.emplace_back(mutant, num_mutants);

  PrivateRepGame prg(pop, params["seed"].get<uint64_t>());
  prg.Update(params.at("t_init"), params.at("q"), params.at("mu_impl"),
             params.at("mu_assess"), params.at("tau"),
             false);
  prg.ResetCounts();

  prg.Update(params.at("t_measure"), params.at("q"), params.at("mu_impl"),
             params.at("mu_assess"), params.at("tau"), false);
  auto pcs = prg.NormCooperationLevels();
  CooperationProbs output{ .pc_res_res = pcs[0][0], .pc_res_mut = pcs[0][1], .pc_mut_res = pcs[1][0], .pc_mut_mut = pcs[1][1] };

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
    std::cerr << "Usage: " << argv[0] << " [-j param.json] norm1 norm2 num_samples" << std::endl;
    std::cerr << "Default parameters:" << std::endl;
    std::cerr << "  " << default_params.dump(2) << std::endl;
  };


  if (args.size() != 3) {
    std::cerr << "[Error] wrong input format" << std::endl;
    show_usage();
    return 1;
  }
  Norm resident = Norm::ParseNormString(args.front(), false);
  args.pop();
  Norm mutant = Norm::ParseNormString(args.front(), false);
  args.pop();
  size_t num_samples = std::strtoull(args.front().c_str(), nullptr, 10);

  auto start = std::chrono::high_resolution_clock::now();

  size_t N = params.at("N").get<size_t>();

  for (size_t num_mut = 2; num_mut < N; num_mut+=2) {
    std::vector<double> v_pc_res_res(num_samples, 0.0), v_pc_res_mut(num_samples, 0.0), v_pc_mut_res(num_samples, 0.0), v_pc_mut_mut(num_samples, 0.0);
    #pragma omp parallel for
    for (size_t i = 0; i < num_samples; i++) {
      nlohmann::json params_copy = params;  // copy parameter
      params_copy["seed"] = params["seed"].get<uint64_t>() + i;
      auto output = RunPolymorphicSimulation(resident, mutant, num_mut, params_copy, false);
      v_pc_res_res[i] = output.pc_res_res;
      v_pc_res_mut[i] = output.pc_res_mut;
      v_pc_mut_res[i] = output.pc_mut_res;
      v_pc_mut_mut[i] = output.pc_mut_mut;
    }

    // average and error
    auto [pc_res_res_avg, pc_res_res_err] = AverageAndError(v_pc_res_res);
    auto [pc_res_mut_avg, pc_res_mut_err] = AverageAndError(v_pc_res_mut);
    auto [pc_mut_res_avg, pc_mut_res_err] = AverageAndError(v_pc_mut_res);
    auto [pc_mut_mut_avg, pc_mut_mut_err] = AverageAndError(v_pc_mut_mut);
    std::cout << num_mut << " "
      << pc_res_res_avg << " " << pc_res_res_err << " "
      << pc_res_mut_avg << " " << pc_res_mut_err << " "
      << pc_mut_res_avg << " " << pc_mut_res_err << " "
      << pc_mut_mut_avg << " " << pc_mut_mut_err << std::endl;
  }

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cerr << "Elapsed time: " << elapsed.count() << " s\n";

  return 0;
}