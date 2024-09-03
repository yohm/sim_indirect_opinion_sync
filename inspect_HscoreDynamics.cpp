#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <queue>
#include <string>
#include <cerrno>
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
  const nlohmann::json default_params = { {"mu_impl", 0.0}, {"mu_assess", 0.02}, {"dt", 0.01}, {"interval", 100}, {"t_max", 10000} };
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

  if (args.size() != 1 && args.size() != 3) {
    std::cerr << "Usage: " << argv[0] << " resident_norm [mutant_norm nu]" << std::endl;
    return 1;
  }

  double dt = params["dt"];
  size_t interval = params["interval"], t_max = params["t_max"];
  double mu_e = params["mu_impl"], mu_a = params["mu_assess"];

  if (args.size() == 1) {
    Norm res = Norm::ParseNormString(args.front());
    res = res.RescaleWithError(mu_e, mu_a);
    auto hvec = SolitaryObsGame::CalculateHDynamics(res, 0.5, interval, dt, t_max);
    for (size_t i = 0; i < hvec.size(); ++i) {
      std::cout << i * interval * dt << " " << hvec.at(i) << std::endl;
    }
    return 0;
  }

  Norm res = Norm::ParseNormString(args[0]);
  Norm mut = Norm::ParseNormString(args[1]);

  char* endptr = nullptr;
  errno = 0;
  double nu = std::strtod(args[2].c_str(), &endptr);
  if (endptr == argv[3] || errno == ERANGE) {
    std::cerr << "Invalid nu" << std::endl;
    return 1;
  }

  res = res.RescaleWithError(mu_e, mu_a);
  mut = mut.RescaleWithError(mu_e, mu_a);

  auto start = std::chrono::high_resolution_clock::now();

  auto ht = SolitaryObsGame::CalculateHDynamicsPolymorphic(res, mut, nu, {0.0,0.0,0.0,0.0}, interval, dt, t_max);
  std::cerr << nu << " " << ht.size() << std::endl;
  for (size_t i = 0; i < ht.size(); ++i) {
    auto [h_rr, h_rm, h_mr, h_mm] = ht.at(i);
    std::cout << i * interval * dt << " " << h_rr << " " << h_rm << " " << h_mr << " " << h_mm << std::endl;
  }

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cerr << "Elapsed time: " << elapsed.count() << " s\n";

  return 0;
}