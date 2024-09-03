#ifndef NORM_TERNARY_HPP
#define NORM_TERNARY_HPP

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

enum class ReputationT {
  B = 0,   // bad
  N = 1,   // neutral
  G = 2    // good
};

char A2C(Action a) {
  if (a == Action::D) { return 'd'; }
  else { return 'c'; }
}
Action C2A(char c) {
  if (c == 'd') { return Action::D; }
  else if (c == 'c') { return Action::C; }
  else { throw std::runtime_error("invalid character for action"); }
}
char R2C(ReputationT r) {
  if (r == ReputationT::B) { return 'B'; }
  else if (r == ReputationT::N) { return 'N'; }
  else { return 'G'; }
}
ReputationT C2R(char c) {
  if (c == 'B') { return ReputationT::B; }
  else if (c == 'N') { return ReputationT::N; }
  else if (c == 'G') { return ReputationT::G; }
  else { throw std::runtime_error("invalid character for action"); }
}
std::ostream &operator<<(std::ostream &os, const Action &act) {
  os << A2C(act);
  return os;
}
std::ostream &operator<<(std::ostream &os, const ReputationT &rep) {
  os << R2C(rep);
  return os;
}

class ActionRuleT {
  public:
  explicit ActionRuleT(double pc_G, double pc_N, double pc_B) : pc_G(pc_G), pc_N(pc_N), pc_B(pc_B) {};
  // pc_B: cooperation probability against Bad recipient, pc_G: cooperation probability against Good recipient

  double pc_G, pc_N, pc_B;

  ActionRuleT Clone() const { return ActionRuleT(pc_G, pc_N, pc_B); }

  double CProb(ReputationT rep_recip) const {
    if (rep_recip == ReputationT::B) { return pc_B; }
    else if (rep_recip == ReputationT::N) { return pc_N; }
    else { return pc_G; }
  }
  void SetCProb(ReputationT rep_recip, double c_prob) {
    if (rep_recip == ReputationT::G) { pc_G = c_prob; }
    else if (rep_recip == ReputationT::B) { pc_B = c_prob; }
    else { pc_N = c_prob; }
  }

  std::string Inspect() const {
    std::stringstream ss;
    std::string name;
    if (IsDeterministic()) {
      if (pc_B == 1.0 && pc_N == 1.0 && pc_G == 1.0) { name = "ALLC"; }
      else if (pc_B == 0.0 && pc_N == 0.0 && pc_G == 0.0) { name = "ALLD"; }
    }
    ss << "ActionRuleT: " << name << std::endl;
    ss << "  (->G): " << pc_G << std::endl;
    ss << "  (->N): " << pc_N << std::endl;
    ss << "  (->B): " << pc_B << std::endl;
    return ss.str();
  }

  ActionRuleT RescaleWithError(double mu_e) const {
    return ActionRuleT{pc_G * (1.0 - mu_e), pc_N * (1.0 - mu_e), pc_B * (1.0 - mu_e)};
  }

  bool IsDeterministic() const {
    return ((pc_B == 0.0 || pc_B == 1.0) && (pc_N == 0.0 || pc_N == 1.0) && (pc_G == 0.0 || pc_G == 1.0));
  }

  int ID() const {
    if (!IsDeterministic()) { return -1; }
    int id = 0;
    if (pc_G == 1.0) { id += 4; }
    if (pc_N == 1.0) { id += 2; }
    if (pc_B == 1.0) { id += 1; }
    return id;
  }

  static ActionRuleT MakeDeterministicRule(int id) {
    if (id < 0 || id >= 8) {
      throw std::runtime_error("ActionRuleT: id must be between 0 and 7");
    }
    double pc_G = (id & 4) ? 1.0 : 0.0;
    double pc_N = (id & 2) ? 1.0 : 0.0;
    double pc_B = (id & 1) ? 1.0 : 0.0;
    return ActionRuleT{pc_G, pc_N, pc_B};
  }

  static ActionRuleT ALLC() {
    return ActionRuleT{1.0, 1.0, 1.0};
  }
  static ActionRuleT ALLD() {
    return ActionRuleT{0.0, 0.0, 0.0};
  }
};

bool operator==(const ActionRuleT& t1, const ActionRuleT& t2) {
  return t1.pc_B == t2.pc_B && t1.pc_N == t2.pc_N && t1.pc_G == t2.pc_G;
}
bool operator!=(const ActionRuleT& t1, const ActionRuleT& t2) { return !(t1 == t2); }

class AssessmentRuleT {
  public:
  explicit AssessmentRuleT(ReputationT r_gc, ReputationT r_gd, ReputationT r_nc, ReputationT r_nd, ReputationT r_bc, ReputationT r_bd)
      : r_gc(r_gc), r_gd(r_gd), r_nc(r_nc), r_nd(r_nd), r_bc(r_bc), r_bd(r_bd) {};
  ReputationT r_gc, r_gd, r_nc, r_nd, r_bc, r_bd;
  // We only consider deterministic assessment rules

  AssessmentRuleT Clone() const { return AssessmentRuleT{r_gc, r_gd, r_nc, r_nd, r_bc, r_bd}; }
  ReputationT Assess(ReputationT rep_recip, Action act) const {
    if (rep_recip == ReputationT::B) {
      if (act == Action::C) { return r_bc; }
      else { return r_bd; }
    }
    else if (rep_recip == ReputationT::N) {
      if (act == Action::C) { return r_nc; }
      else { return r_nd; }
    }
    else {
      if (act == Action::C) { return r_gc; }
      else { return r_gd; }
    }
  }
  void SetAssessment(ReputationT rep_recip, Action act, ReputationT new_rep) {
    if (rep_recip == ReputationT::B) {
      if (act == Action::C) { r_bc = new_rep; }
      else { r_bd = new_rep; }
    }
    else if (rep_recip == ReputationT::N) {
      if (act == Action::C) { r_nc = new_rep; }
      else { r_nd = new_rep; }
    }
    else {
      if (act == Action::C) { r_gc = new_rep; }
      else { r_gd = new_rep; }
    }
  }

  std::string Inspect() const {
    std::stringstream ss;
    ss << "AssessmentRuleT: " << std::endl;
    ss << "  (G->C): " << r_gc << std::endl;
    ss << "  (G->D): " << r_gd << std::endl;
    ss << "  (N->C): " << r_nc << std::endl;
    ss << "  (N->D): " << r_nd << std::endl;
    ss << "  (B->C): " << r_bc << std::endl;
    ss << "  (B->D): " << r_bd << std::endl;
    return ss.str();
  }

  int ID() const {
    int id = 0;
    for (ReputationT r: std::array<ReputationT,6>{r_gc, r_gd, r_nc, r_nd, r_bc, r_bd}) {
      id *= 3;
      int r_int = static_cast<int>(r);
      id += r_int;
    }
    return id;
  }

  static AssessmentRuleT MakeDeterministicRule(int id) {
    if (id < 0 || id >= 729) {
      throw std::runtime_error("AssessmentRuleT: id must be between 0 and 729");
    }
    int r_bd = id % 3;
    id /= 3;
    int r_bc = id % 3;
    id /= 3;
    int r_nd = id % 3;
    id /= 3;
    int r_nc = id % 3;
    id /= 3;
    int r_gd = id % 3;
    id /= 3;
    int r_gc = id % 3;
    return AssessmentRuleT{static_cast<ReputationT>(r_gc), static_cast<ReputationT>(r_gd),
                           static_cast<ReputationT>(r_nc), static_cast<ReputationT>(r_nd),
                           static_cast<ReputationT>(r_bc), static_cast<ReputationT>(r_bd)};
  }

  static AssessmentRuleT AllGood() {
    return AssessmentRuleT::MakeDeterministicRule(728);
  }
  static AssessmentRuleT AllBad() {
    return AssessmentRuleT::MakeDeterministicRule(0);
  }
};

bool operator==(const AssessmentRuleT& t1, const AssessmentRuleT& t2) {
  return t1.r_gc == t2.r_gc && t1.r_gd == t2.r_gd && t1.r_nc == t2.r_nc &&
         t1.r_nd == t2.r_nd && t1.r_bc == t2.r_bc && t1.r_bd == t2.r_bd;
}
bool operator!=(const AssessmentRuleT& t1, const AssessmentRuleT& t2) {
  return !(t1 == t2);
}

// NormT is a set of AssessmentRuleT & ActionRuleT
class NormT {
public:
  NormT(const AssessmentRuleT &R, const ActionRuleT &P)
      : R(R), P(P) {};
  NormT(const NormT &rhs) : R(rhs.R), P(rhs.P) {};
  NormT(ReputationT r_gc, ReputationT r_gd, ReputationT r_nc, ReputationT r_nd, ReputationT r_bc, ReputationT r_bd,
        double pc_G, double pc_N, double pc_B)
      : R(r_gc, r_gd, r_nc, r_nd, r_bc, r_bd), P(pc_G, pc_N, pc_B) {};
  AssessmentRuleT R;  // assessment rules to assess donor and recipient
  ActionRuleT P;  // action rule
  std::string Inspect() const {
    std::stringstream ss;
    ss << "NormT: " << std::endl;
    if (IsDeterministic()) {
      ss << "  ID: " << std::bitset<6>(ID()) << std::endl;
    }
    for (ReputationT r: std::array<ReputationT,3>{ReputationT::G, ReputationT::N, ReputationT::B}) {
      double c_prob = P.CProb(r);
      ReputationT rep_c = R.Assess(r, Action::C);
      ReputationT rep_d = R.Assess(r, Action::D);
      ss << std::setprecision(3) << std::fixed;
      ss << "( ->" << r << "): \n"
         << "  P:" << c_prob << "\n"
         << "  R (c:" << rep_c << ",d:" << rep_d << ")\n";

    }
    return ss.str();
  }
  double CProb(ReputationT recipient) const { return P.CProb(recipient); }
  ReputationT Assess(ReputationT recipient, Action act) const {
    return R.Assess(recipient, act);
  }
  bool IsDeterministic() const {
    return P.IsDeterministic();
  }
  int ID() const {
    if (!IsDeterministic()) {
      return -1;
    }
    int id = 0;
    id += R.ID() << 3;
    id += P.ID();
    return id;
  }
  static NormT ConstructFromID(int id) {
    if (id < 0 || id >= 729 * 8) {
      throw std::runtime_error("NormT: id must be between 0 and 5831");
    }
    int R_id = id >> 3;
    int P_id = id & 0b111;
    return NormT{AssessmentRuleT::MakeDeterministicRule(R_id),
                 ActionRuleT::MakeDeterministicRule(P_id)};
  }
  static NormT AllC() {
    return NormT{AssessmentRuleT::MakeDeterministicRule(728),
                 ActionRuleT{1, 1, 1}};
  }
  static NormT AllD() {
    return NormT{AssessmentRuleT::MakeDeterministicRule(0),
                 ActionRuleT{0, 0, 0}};
  }
  static NormT Random(double p = 0.5) {
    return NormT{ AssessmentRuleT::MakeDeterministicRule(0),
                  ActionRuleT{p, p, p} };
  }

  std::string GetName() const {
    std::ostringstream oss;
    oss << R2C(R.r_gc) << R2C(R.r_gd) << R2C(R.r_nc) << R2C(R.r_nd) << R2C(R.r_bc) << R2C(R.r_bd) << " "
        << std::setprecision(3) << std::fixed << P.pc_G << " " << P.pc_N << " " << P.pc_B;
    return oss.str();
  }

  static NormT ParseNormString(const std::string& str) {
    std::regex re_d(R"(^\d+$)"); // regex for digits
    std::regex re_x(R"(^0x[0-9a-fA-F]+$)");  // regex for digits in hexadecimal
    // regular expression for 9 floating point numbers separated by space
    std::regex re_a(R"(^[GNB]{6} (?:\d+|\d+\.\d+) (?:\d+|\d+\.\d+) (?:\d+|\d+\.\d+)$)");
    NormT norm = NormT::AllC();
    if (std::regex_match(str, re_d)) {
      int id = std::stoi(str);
      norm = NormT::ConstructFromID(id);
    }
    else if (std::regex_match(str, re_x)) {
      int id = std::stoi(str, nullptr, 16);
      norm = NormT::ConstructFromID(id);
    }
    else if (std::regex_match(str, re_a)) {
      std::istringstream iss(str);
      std::string rep_str;
      iss >> rep_str;
      AssessmentRuleT R{C2R(rep_str[0]), C2R(rep_str[1]), C2R(rep_str[2]), C2R(rep_str[3]), C2R(rep_str[4]), C2R(rep_str[5])};
      double p_G, p_N, p_B;
      iss >> p_G >> p_N >> p_B;
      // p_G, p_N, p_B must be in [0, 1]
      if (p_G < 0.0 || p_G > 1.0 || p_N < 0.0 || p_N > 1.0 || p_B < 0.0 || p_B > 1.0) {
        std::cerr << "[Error] invalid probability: " << str << std::endl;
        throw std::runtime_error("invalid probability");
      }
      ActionRuleT P{p_G, p_N, p_B};
      norm = NormT(R, P);
    }
    else {
      std::cerr << "[Error] invalid norm string: " << str << std::endl;
      throw std::runtime_error("invalid norm string");
    }
    return norm;
  }
};

bool operator==(const NormT& n1, const NormT& n2) {
  return n1.P == n2.P && n1.R == n2.R;
}
bool operator!=(const NormT& n1, const NormT& n2) {
  return !(n1 == n2);
}

#endif
