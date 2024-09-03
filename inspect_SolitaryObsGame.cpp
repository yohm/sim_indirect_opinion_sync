#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <queue>
#include <string>
#include <nlohmann/json.hpp>
#include "SolitaryObsGame.hpp"


int main(int argc, char *argv[]) {

  std::vector<std::string> args;
  nlohmann::json params = nlohmann::json::object();
  // -j param.json : set parameters by json file or json string
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "-j" && i + 1 < argc) {
      std::ifstream fin(argv[++i]);
      // check if file exists
      if (fin) {
        fin >> params;
        fin.close();
      }
      else {
        std::istringstream iss(argv[i]);
        iss >> params;
      }
    }
    else {
      args.emplace_back(argv[i]);
    }
  }

  // set default parameters
  const nlohmann::json default_params = { {"mu_impl", 0.0}, {"mu_assess", 0.02} };
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

  if (args.size() != 1 && args.size() != 2 && args.size() != 3) {
    std::cerr << "Usage: " << argv[0] << " [-j param.json] resident_norm_id [mutant_norm_id] [nu]" << std::endl;
    return 1;
  }

  double mu_a = params["mu_assess"];
  double mu_e = params["mu_impl"];

  if (args.size() == 1) {
    Norm res = Norm::ParseNormString(args[0]);
    res = res.RescaleWithError(mu_e, mu_a);
    double h = SolitaryObsGame::HscoreMonomorphic(res);
    double pc = SolitaryObsGame::SelfCooperationLevel(res, h);
    std::cout << "h: " << h << " pc: " << pc << std::endl;
    return 0;
  }

  Norm res = Norm::ParseNormString(args[0]);
  Norm mut = Norm::ParseNormString(args[1]);

  res = res.RescaleWithError(mu_e, mu_a);
  mut = mut.RescaleWithError(mu_e, mu_a);

  if (args.size() == 3) {
    // parse nu
    char* endptr = nullptr;
    errno = 0;
    double nu = std::strtod(args[2].c_str(), &endptr);
    if (endptr == argv[3] || errno == ERANGE) {
      std::cerr << "Invalid nu" << std::endl;
      return 1;
    }
    auto [pc_rr,pc_rm,pc_mr,pc_mm] = SolitaryObsGame::CooperationLevelsPolymorphic(res, mut, nu);
    std::cout << "pc_rr: " << pc_rr << " pc_rm: " << pc_rm << std::endl << "pc_mr: " << pc_mr << " pc_mm: " << pc_mm << std::endl;
    return 0;
  }

  double benefit = 5.0;

  auto start = std::chrono::high_resolution_clock::now();

  for (int n = 0; n < 50; n++) {
    double nu = n * 0.02;
    auto [pi_res,pi_mut] = SolitaryObsGame::Payoffs(res, mut, nu, benefit);
    std::cout << nu << " " << pi_res << " " << pi_mut << std::endl;
  }

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cerr << "Elapsed time: " << elapsed.count() << " s\n";

  return 0;
}