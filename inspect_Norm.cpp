#include <iostream>
#include <bitset>
#include <vector>
#include <set>
#include "Norm.hpp"

constexpr Reputation B = Reputation::B, G = Reputation::G;
constexpr Action C = Action::C, D = Action::D;

int main(int argc, char** argv) {

  std::vector<std::string> args;
  bool swap_gb = false;
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "-s") {
      swap_gb = true;
    }
    else {
      args.emplace_back(argv[i]);
    }
  }

  if (args.size() == 1) {
    Norm n = Norm::ParseNormString(args[0], swap_gb);
    std::cout << n.Inspect();
  }
  else if (args.size() >= 2) {
    Norm n = Norm::ParseNormString(args[0], swap_gb);
    // loop over the other norms
  }
  else {   // no arguments
    std::cerr << "Usage: " << argv[0] << " [options] norm [other norms]" << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << "  -s         : swap G and B" << std::endl;
    std::cerr << "Norm format: [Norm name] or [ID] or [R_BBd, R_BBc, R_BGd, R_BGc, R_GBd, R_GBc, R_GGd, R_GGc, P_BB, P_BG, P_GB, P_GG]" << std::endl;
    return 0;
  }

  return 0;
}