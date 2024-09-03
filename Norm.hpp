#ifndef NORM_HPP
#define NORM_HPP

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <array>
#include <map>
#include <algorithm>
#include <cstdint>
#include <regex>


enum class Action {
  D = 0,  // defect
  C = 1   // cooperate
};

Action FlipAction(Action a) {
  if (a == Action::C) { return Action::D; }
  else { return Action::C; }
}

enum class Reputation {
  B = 0,   // bad
  G = 1    // good
};

Reputation FlipReputation(Reputation r) {
  if (r == Reputation::G) { return Reputation::B; }
  else { return Reputation::G; }
}

char A2C(Action a) {
  if (a == Action::D) { return 'd'; }
  else { return 'c'; }
}
Action C2A(char c) {
  if (c == 'd') { return Action::D; }
  else if (c == 'c') { return Action::C; }
  else { throw std::runtime_error("invalid character for action"); }
}
char R2C(Reputation r) {
  if (r == Reputation::B) { return 'B'; }
  else { return 'G'; }
}
Reputation C2R(char c) {
  if (c == 'B') { return Reputation::B; }
  else if (c == 'G') { return Reputation::G; }
  else { throw std::runtime_error("invalid character for action"); }
}
std::ostream &operator<<(std::ostream &os, const Action &act) {
  os << A2C(act);
  return os;
}
std::ostream &operator<<(std::ostream &os, const Reputation &rep) {
  os << R2C(rep);
  return os;
}

class ActionRule {
  public:
  explicit ActionRule(double pc_G, double pc_B) : pc_G(pc_G), pc_B(pc_B) {};
  // pc_B: cooperation probability against Bad recipient, pc_G: cooperation probability against Good recipient

  double pc_G, pc_B;

  ActionRule Clone() const { return ActionRule(pc_G, pc_B); }

  double CProb(Reputation rep_recip) const {
    return (rep_recip == Reputation::B) ? pc_B : pc_G;
  }
  void SetCProb(Reputation rep_recip, double c_prob) {
    if (rep_recip == Reputation::B) { pc_B = c_prob; }
    else { pc_G = c_prob; }
  }
  ActionRule SwapGB() const { // returns a new action rule with G and B swapped
    return ActionRule{pc_B, pc_G};
  }

  std::string Inspect() const {
    std::stringstream ss;
    std::string name;
    if (IsDeterministic()) {
      if (pc_B == 1.0 && pc_G == 1.0) { name = "ALLC"; }
      else if (pc_B == 0.0 && pc_G == 0.0) { name = "ALLD"; }
      else if (pc_B == 0.0 && pc_G == 1.0) { name = "DISC"; }
      else if (pc_B == 1.0 && pc_G == 0.0) { name = "ADISC"; }
    }
    ss << "ActionRule: " << name << std::endl;
    ss << "  (->G): " << pc_G << std::endl;
    ss << "  (->B): " << pc_B << std::endl;
    return ss.str();
  }

  ActionRule RescaleWithError(double mu_e) const {
    return ActionRule{pc_G * (1.0 - mu_e), pc_B * (1.0 - mu_e)};
  }

  bool IsDeterministic() const {
    return ((pc_B == 0.0 || pc_B == 1.0) && (pc_G == 0.0 || pc_G == 1.0));
  }

  int ID() const {
    if (!IsDeterministic()) { return -1; }
    int id = 0;
    if (pc_G == 1.0) { id += 2; }
    if (pc_B == 1.0) { id += 1; }
    return id;
  }

  static ActionRule MakeDeterministicRule(int id) {
    if (id == 0) {
      return ActionRule{0.0, 0.0};
    }
    else if (id == 1) {
      return ActionRule{0.0, 1.0};
    }
    else if (id == 2) {
      return ActionRule{1.0, 0.0};
    }
    else if (id == 3) {
      return ActionRule{1.0, 1.0};
    }
    else {
      throw std::runtime_error("ActionRuleDet: id must be between 0 and 3");
    }
  }

  static ActionRule DISC() {
    return ActionRule{1.0, 0.0};
  }
  static ActionRule ADISC() {  // anti-discriminator
    return ActionRule{0.0, 1.0};
  }
  static ActionRule ALLC() {
    return ActionRule{1.0, 1.0};
  }
  static ActionRule ALLD() {
    return ActionRule{0.0, 0.0};
  }
};

bool operator==(const ActionRule& t1, const ActionRule& t2) {
  return t1.pc_B == t2.pc_B && t1.pc_G == t2.pc_G;
}
bool operator!=(const ActionRule& t1, const ActionRule& t2) { return !(t1 == t2); }

class AssessmentRule {
  public:
  explicit AssessmentRule(double r_gc, double r_gd, double r_bc, double r_bd)
    : r_gc(r_gc), r_gd(r_gd), r_bc(r_bc), r_bd(r_bd) {};

  double r_gc, r_gd, r_bc, r_bd;

  AssessmentRule Clone() const { return AssessmentRule{r_gc, r_gd, r_bc, r_bd}; }
  AssessmentRule SwapGB() const {
    return AssessmentRule{ 1.0-r_bc, 1.0-r_bd, 1.0-r_gc, 1.0-r_gd };
  }
  double GProb(Reputation rep_recip, Action act) const {
    if (rep_recip == Reputation::B) {
      return (act == Action::C) ? r_bc : r_bd;
    }
    else {
      return (act == Action::C) ? r_gc : r_gd;
    }
  }
  void SetGProb(Reputation rep_recip, Action act, double g_prob) {
    if (rep_recip == Reputation::B) {
      if (act == Action::C) { r_bc = g_prob; }
      else { r_bd = g_prob; }
    }
    else {
      if (act == Action::C) { r_gc = g_prob; }
      else { r_gd = g_prob; }
    }
  }

  std::string Inspect() const {
    std::stringstream ss;
    ss << "AssessmentRule: " << std::endl;
    ss << "  (G->C): " << r_gc << std::endl;
    ss << "  (G->D): " << r_gd << std::endl;
    ss << "  (B->C): " << r_bc << std::endl;
    ss << "  (B->D): " << r_bd << std::endl;
    return ss.str();
  }

  AssessmentRule RescaleWithError(double mu_a) const {
    return AssessmentRule{ r_gc * (1.0 - 2.0 * mu_a) + mu_a,
                           r_gd * (1.0 - 2.0 * mu_a) + mu_a,
                           r_bc * (1.0 - 2.0 * mu_a) + mu_a,
                           r_bd * (1.0 - 2.0 * mu_a) + mu_a };
  }

  bool IsDeterministic() const {
    return ((r_gc == 0.0 || r_gc == 1.0) && (r_gd == 0.0 || r_gd == 1.0) &&
            (r_bc == 0.0 || r_bc == 1.0) && (r_bd == 0.0 || r_bd == 1.0));
  }

  int ID() const {
    if (!IsDeterministic()) { return -1; }
    int id = 0;
    if (r_gc == 1.0) { id += 8; }
    if (r_gd == 1.0) { id += 4; }
    if (r_bc == 1.0) { id += 2; }
    if (r_bd == 1.0) { id += 1; }
    return id;
  }

  static AssessmentRule MakeDeterministicRule(int id) {
    if (id < 0 || id >= 16) {
      throw std::runtime_error("AssessmentRule: id must be between 0 and 15");
    }
    double r_gc = (id & 8) ? 1.0 : 0.0;
    double r_gd = (id & 4) ? 1.0 : 0.0;
    double r_bc = (id & 2) ? 1.0 : 0.0;
    double r_bd = (id & 1) ? 1.0 : 0.0;
    return AssessmentRule{r_gc, r_gd, r_bc, r_bd};
  }

  static AssessmentRule AllGood() {
    return AssessmentRule::MakeDeterministicRule(0b1111);
  }
  static AssessmentRule AllBad() {
    return AssessmentRule::MakeDeterministicRule(0b0000);
  }
  static AssessmentRule ImageScoring() {
    return AssessmentRule{1.0, 0.0, 1.0, 0.0};
  }
};

bool operator==(const AssessmentRule& t1, const AssessmentRule& t2) {
  return t1.r_gc == t2.r_gc && t1.r_gd == t2.r_gd && t1.r_bc == t2.r_bc && t1.r_bd == t2.r_bd;
}
bool operator!=(const AssessmentRule& t1, const AssessmentRule& t2) {
  return !(t1 == t2);
}

// Norm is a set of AssessmentRule & ActionRule
class Norm {
public:
  Norm(const AssessmentRule &R, const ActionRule &P)
      : R(R), P(P) {};
  Norm(const Norm &rhs) : R(rhs.R), P(rhs.P) {};
  Norm(double r_gc, double r_gd, double r_bc, double r_bd, double pc_G, double pc_B)
      : R(r_gc, r_gd, r_bc, r_bd), P(pc_G, pc_B) {};
  AssessmentRule R;  // assessment rules to assess donor and recipient
  ActionRule P;  // action rule
  std::string Inspect() const {
    std::stringstream ss;
    ss << "Norm: " << std::endl;
    if (IsDeterministic()) {
      ss << "  ID: " << std::bitset<6>(ID()) << std::endl;
    }
    for (Reputation r: std::array<Reputation,2>{Reputation::G, Reputation::B}) {
      double c_prob = P.CProb(r);
      double g_prob_c = R.GProb(r, Action::C);
      double g_prob_d = R.GProb(r, Action::D);
      ss << std::setprecision(3) << std::fixed;
      ss << "( ->" << r << "): \n"
         << "  P:" << c_prob << "\n"
         << "  R (c:" << g_prob_c << ",d:" << g_prob_d << ")\n";

    }
    return ss.str();
  }
  double CProb(Reputation recipient) const { return P.CProb(recipient); }
  double Gprob(Reputation recipient, Action act) const {
    return R.GProb(recipient, act);
  }
  Norm SwapGB() const {
    return Norm{R.SwapGB(), P.SwapGB()};
  }
  bool IsDeterministic() const {
    return P.IsDeterministic() && R.IsDeterministic();
  }
  Norm RescaleWithError(double mu_e, double mu_a) const {
    return Norm{R.RescaleWithError(mu_a), P.RescaleWithError(mu_e)};
  }
  int ID() const {
    if (!IsDeterministic()) {
      return -1;
    }
    int id = 0;
    id += R.ID() << 2;
    id += P.ID();
    return id;
  }
  static Norm ConstructFromID(int id) {
    if (id < 0 || id >= (1 << 6)) {
      throw std::runtime_error("Norm: id must be between 0 and 63");
    }
    int R_id = (id >> 2) & 0b1111;
    int P_id = id & 0b11;
    return Norm{AssessmentRule::MakeDeterministicRule(R_id),
                ActionRule::MakeDeterministicRule(P_id)};
  }
  static Norm AllC() {
    return Norm{AssessmentRule{1, 1, 1, 1},
                ActionRule{1, 1}};
  }
  static Norm AllD() {
    return Norm{AssessmentRule{0, 0, 0, 0},
                ActionRule{0, 0}};
  }
  static Norm AllG() {
    // Always assess G, but action rule is Discriminator
    // it may defect under assessment error
    return Norm{AssessmentRule::AllGood(),
                ActionRule::DISC()};
  }
  static Norm AllB() {
    // Always assess B, but action rule is Discriminator
    // it may defect under assessment error
    return Norm{AssessmentRule::AllBad(),
                ActionRule::DISC()};
  }
  static Norm ImageScoring() {
    return Norm{AssessmentRule::ImageScoring(),
                ActionRule::DISC()};
  }
  static Norm Random(double p = 0.5) {
    return Norm{ AssessmentRule{0, 0, 0, 0},
                 ActionRule{p, p} };
  }
  static Norm L3() {
    return Norm{AssessmentRule{1, 0, 1, 1},
                ActionRule{1, 0}};
  }
  static Norm L6() {
    return Norm{AssessmentRule{1, 0, 0, 1},
                ActionRule{1, 0}};
  }
  static Norm GenerousScoring(double benefit = 5.0) {
    double q = 1.0 / benefit;
    AssessmentRule gsco(1.0, q, 1.0, q);
    return Norm{gsco, ActionRule::DISC()};
  }
  bool IsGenerousScoring() const {
    if( P != ActionRule::DISC() ) {
      return false;
    }
    if (R.r_gc != 1.0 || R.r_bc != 1.0) {
      return false;
    }
    if (R.r_gd != R.r_bd) {
      return false;
    }
    return true;
  }
  static std::vector<Norm> DeterministicNorms() {
    std::vector<Norm> norms;
    for (int i = 0; i < 64; i++) {
      const Norm n = Norm::ConstructFromID(i);
      const Norm ns = n.SwapGB();
      int nid = (n.P.ID() << 4) + n.R.ID();
      int nsid = (ns.P.ID() << 4) + ns.R.ID();
      if (nid >= nsid) {
        norms.push_back(n);
      }
    }
    return norms;
  }
  static const std::vector<std::pair<int, std::string> > NormNames;
  std::string GetName() const {
    if (IsDeterministic()) {
      int id = ID();
      for (auto &p : NormNames) {
        if (p.first == id) {
          return p.second;
        }
      }
    }
    if (IsGenerousScoring()) {
      double benefit = 1.0 / R.r_gd;
      std::ostringstream oss;
      oss << "GSCO-" << std::fixed << std::setprecision(1) << benefit;
      return oss.str();
    }
    return "";
  }
  static Norm ConstructFromName(const std::string &name) {
    auto found = std::find_if(Norm::NormNames.begin(), Norm::NormNames.end(), [&](auto &s) {
      return s.second == name;
    });

    std::regex re_gsco(R"(^GSCO-(\d+(?:\.\d+)?)$)");
    std::regex re_second(R"(^SECOND-(\d+)$)"); // regex for SECOND-[id] (id: 0-15
    std::regex re_random(R"(^RANDOM-(0(\.\d+)?|1(\.0+)?)$)");  // RANDOM-0.5
    std::smatch matches;

    if (found != Norm::NormNames.end()) {
      return Norm::ConstructFromID(found->first);
    }
    else if (std::regex_match(name, matches, re_gsco)) {
      double benefit = std::stod(matches[1]);
      return Norm::GenerousScoring(benefit);
    }
    else if (std::regex_match(name, matches, re_second)) {
      int id = std::stoi(matches[1]);
      if (id < 0 || id >= 16) {
        throw std::runtime_error("SECOND: id must be between 0 and 15");
      }
      return {AssessmentRule::MakeDeterministicRule(id), ActionRule::DISC()};
    }
    else if (std::regex_match(name, matches, re_random)) {
      double p = std::stod(matches[1]);
      return Norm::Random(p);
    }
    else if (name == "RANDOM") {
      return Norm::Random(0.5);
    }
    else {
      std::cerr << "Unknown norm: " << name << std::endl;
      std::cerr << "Available norms are: " << std::endl;
      for (auto &n : Norm::NormNames) {
        std::cerr << "\t" << n.second << std::endl;
      }
      std::cerr << "\t" << "GSCO-benefit" << std::endl;
      std::cerr << "\t" << "SECOND-[id] (id: 0-15)" << std::endl;
      std::cerr << "\t" << "RANDOM-[p] (p: a real number in [0, 1])" << std::endl;
      throw std::runtime_error("Unknown norm");
    }
  }

  static Norm ParseNormString(const std::string& str, bool swap_gb = false) {
    std::regex re_d(R"(\d+)"); // regex for digits
    std::regex re_x(R"(^0x[0-9a-fA-F]+$)");  // regex for digits in hexadecimal
    // regular expression for 6 floating point numbers separated by space
    std::regex re_a(R"(^(0|1)(\.\d+)?( (0|1)(\.\d+)?){5}$)");
    Norm norm = Norm::AllC();
    if (std::regex_match(str, re_d)) {
      int id = std::stoi(str);
      norm = Norm::ConstructFromID(id);
    }
    else if (std::regex_match(str, re_x)) {
      int id = std::stoi(str, nullptr, 16);
      norm = Norm::ConstructFromID(id);
    }
    else if (std::regex_match(str, re_a)) {
      std::istringstream iss(str);
      std::array<double,6> prob{};
      for (int i = 0; i < 6; ++i) {
        iss >> prob[i];
      }
      norm = Norm(prob[0], prob[1], prob[2], prob[3], prob[4], prob[5]);
    }
    else {
      norm = Norm::ConstructFromName(str);
    }
    if (swap_gb) {
      norm = norm.SwapGB();
    }
    return norm;
  }
};

const std::vector<std::pair<int,std::string> > Norm::NormNames = {{
                                                                    {AllC().ID(), "AllC"},
                                                                    {AllD().ID(), "AllD"},
                                                                    {AllG().ID(), "AllG"},
                                                                    {AllB().ID(), "AllB"},
                                                                    {ImageScoring().ID(), "ImageScoring"},
                                                                    {L3().ID(), "L3"},
                                                                    {L6().ID(), "L6"},
                                                                }};
bool operator==(const Norm& n1, const Norm& n2) {
  return n1.P == n2.P && n1.R == n2.R;
}
bool operator!=(const Norm& n1, const Norm& n2) {
  return !(n1 == n2);
}

#endif
