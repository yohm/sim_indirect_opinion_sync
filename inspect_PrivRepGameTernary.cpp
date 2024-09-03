#include <iostream>
#include <fstream>
#include <chrono>
#include <regex>
#include <vector>
#include <queue>
#include <utility>
#include <string>
#include <nlohmann/json.hpp>
#include "PrivRepGameTernary.hpp"


void PrintInitialTimeSeries(PrivateRepGameT& prg, const nlohmann::json& params, std::ostream& out = std::cout) {
  size_t num_prints = 100;
  size_t interval = params.at("t_init").get<size_t>() / num_prints;

  size_t N_strategies = prg.Population().size();  // number of different strategies

  for (size_t t = 0; t < num_prints; t++) {
    prg.Update(interval, params.at("q"), params.at("mu_impl"), params.at("mu_assess"), false);
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
      }
      else {
        std::istringstream iss(argv[i]);
        iss >> params;
      }
    }
    else {
      args.emplace(argv[i]);
    }
  }

  // set default parameters
  const nlohmann::json default_params = { {"t_init", 1e3}, {"t_measure", 1e3}, {"q", 1.0}, {"mu_impl", 0.0}, {"mu_assess", 0.02}, {"seed", 123456789ull} };
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

  auto show_usage = [&argv, default_params] {
    std::cerr << "Usage: " << argv[0] << " [-j param.json] norm1 size1 [norm2 size2 ...]" << std::endl;
    std::cerr << "Default parameters:" << std::endl;
    std::cerr << "  " << default_params.dump(2) << std::endl;
  };

  auto start = std::chrono::high_resolution_clock::now();

  if (args.size() < 2 || args.size() % 2 != 0) {
    std::cerr << "[Error] wrong input format" << std::endl;
    show_usage();
    return 1;
  }
  else {
    PrivateRepGameT::population_t pop;
    while (!args.empty()) {
      std::string norm_str = args.front(); args.pop();
      size_t size = std::stoi(args.front()); args.pop();
      NormT norm = NormT::ParseNormString(norm_str);
      pop.emplace_back(norm, size);
      std::cout << norm.Inspect() << std::endl;
    }

    PrivateRepGameT prg(pop, params["seed"].get<uint64_t>());
    PrintInitialTimeSeries(prg, params, std::cout);

    if (params["t_measure"].get<size_t>() > 0) {
      prg.ResetCounts();
      // calculate goodness for homogeneous population
      bool count_goodness = (pop.size() == 1);
      prg.Update(params["t_measure"].get<size_t>(), params["q"].get<double>(), params["mu_impl"].get<double>(),
                 params["mu_assess"].get<double>(), count_goodness);
      std::cout << "SystemWideCooperationLevel: " << prg.SystemWideCooperationLevel() << std::endl;
      auto c_levels = prg.NormCooperationLevels();
      if (c_levels.size() > 1) {
        std::cout << "NormCooperationLevels:\n";
        for (size_t i = 0; i < c_levels.size(); i++) {
          for (size_t j = 0; j < c_levels[i].size(); j++) {
            std::cout << ' ' << c_levels[i][j];
          }
          std::cout << "\n";
        }
      }
      if (count_goodness) {
        std::cout << "Reputations:\n";
        auto reps = prg.NormAverageReputation();
        for (size_t i = 0; i < reps.size(); i++) {
          for (size_t j = 0; j < reps[i].size(); j++) {
            std::cout << ' ' << reps[i][j][0] << '-' << reps[i][j][1] << '-' << reps[i][j][2];
          }
          std::cout << "\n";
        }
      }
    }

  }

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cerr << "Elapsed time: " << elapsed.count() << " s\n";

  return 0;
}