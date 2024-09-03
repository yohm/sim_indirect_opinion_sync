#ifndef NORM_3RD_HPP
#define NORM_3RD_HPP

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
  ActionRule(const std::array<double,4>& coop_probs) : coop_probs(coop_probs) {};
  // {P(C|B,B), P(C|B,G), P(C|G,B), P(C|G,G)}

  ActionRule(const std::map<std::pair<Reputation,Reputation>,double>& actions) {
    // { (B,B,P(B,B)), (B,G,P(B,G)), (G,B,P(G,B)), (G,G,P(G,G)) }
    if (actions.size() != 4) {
      throw std::runtime_error("unspecified actions");
    }
    coop_probs[0] = actions.at({Reputation::B, Reputation::B});
    coop_probs[1] = actions.at({Reputation::B, Reputation::G});
    coop_probs[2] = actions.at({Reputation::G, Reputation::B});
    coop_probs[3] = actions.at({Reputation::G, Reputation::G});
  }
  ActionRule(const std::map<std::pair<Reputation,Reputation>,Action>& actions) {
    // { (B,B,D), (B,G,C), (G,B,D), (G,G,C) }
    if (actions.size() != 4) {
      throw std::runtime_error("unspecified actions");
    }
    coop_probs[0] = actions.at({Reputation::B, Reputation::B}) == Action::C ? 1.0 : 0.0;
    coop_probs[1] = actions.at({Reputation::B, Reputation::G}) == Action::C ? 1.0 : 0.0;
    coop_probs[2] = actions.at({Reputation::G, Reputation::B}) == Action::C ? 1.0 : 0.0;
    coop_probs[3] = actions.at({Reputation::G, Reputation::G}) == Action::C ? 1.0 : 0.0;
  }

  std::array<double,4> coop_probs; // P(C|B,B), P(C|B,G), P(C|G,B), P(C|G,G)

  ActionRule Clone() const { return ActionRule(coop_probs); }

  double CProb(const Reputation& rep_d, const Reputation& rep_r) const {
    size_t idx = static_cast<size_t>(rep_d) * 2 + static_cast<size_t>(rep_r);
    return coop_probs[idx];
  }
  void SetCProb(const Reputation& rep_d, const Reputation& rep_r, double c_prob) {
    size_t idx = static_cast<size_t>(rep_d) * 2 + static_cast<size_t>(rep_r);
    coop_probs[idx] = c_prob;
  }
  ActionRule SwapGB() const { // returns a new action rule with G and B swapped
    std::array<double,4> new_coop_probs = {coop_probs[3], coop_probs[2], coop_probs[1], coop_probs[0]};
    return ActionRule{new_coop_probs};
  }

  std::string Inspect() const {
    std::stringstream ss;
    ss << "ActionRule: " << ID() << std::endl;
    for (size_t i = 0; i < 4; i++) {
      Reputation rep_d = static_cast<Reputation>(i / 2);
      Reputation rep_r = static_cast<Reputation>(i % 2);
      ss << "(" << rep_d << "->" << rep_r << ") : " << coop_probs[i];
      if (i % 2 == 1) { ss << std::endl; }
      else { ss << "\t"; }
    }
    return ss.str();
  }

  ActionRule RescaleWithError(double mu_e) const {
    std::array<double,4> rescaled = {0.0};
    for (size_t i = 0; i < 4; i++) {
      rescaled[i] = (1.0 - mu_e) * coop_probs[i];
    }
    return ActionRule{rescaled};
  }

  bool IsSecondOrder() const {
    if (coop_probs[0] == coop_probs[2] && coop_probs[1] == coop_probs[3]) { return true; }
    else { return false; }
  }

  bool IsDeterministic() const {
    if (( coop_probs[0] == 0.0 || coop_probs[0] == 1.0 ) &&
        ( coop_probs[1] == 0.0 || coop_probs[1] == 1.0 ) &&
        ( coop_probs[2] == 0.0 || coop_probs[2] == 1.0 ) &&
        ( coop_probs[3] == 0.0 || coop_probs[3] == 1.0 ) ) {
      return true;
    }
    return false;
  }

  int ID() const {
    if (!IsDeterministic()) { return -1; }
    int id = 0;
    for (size_t i = 0; i < 4; i++) {
      if (coop_probs[i] == 1.0) { id += 1 << i; }
    }
    return id;
  }

  static ActionRule MakeDeterministicRule(int id) {
    if (id < 0 || id > 15) {
      throw std::runtime_error("AssessmentRuleDet: id must be between 0 and 15");
    }
    std::array<double,4> c_probs = {0.0};
    for (size_t i = 0; i < c_probs.size(); i++) {
      if ((id >> i) % 2) { c_probs[i] = 1.0; }
      else { c_probs[i] = 0.0; }
    }
    return ActionRule{ c_probs };
  }

  static ActionRule DISC() {
    return ActionRule{ {0.0, 1.0, 0.0, 1.0} };
  }
  static ActionRule ALLC() {
    return ActionRule{ {1.0, 1.0, 1.0, 1.0} };
  }
  static ActionRule ALLD() {
    return ActionRule{ {0.0, 0.0, 0.0, 0.0} };
  }
};

bool operator==(const ActionRule& t1, const ActionRule& t2) {
  return t1.coop_probs[0] == t2.coop_probs[0] &&
         t1.coop_probs[1] == t2.coop_probs[1] &&
         t1.coop_probs[2] == t2.coop_probs[2] &&
         t1.coop_probs[3] == t2.coop_probs[3];
}
bool operator!=(const ActionRule& t1, const ActionRule& t2) { return !(t1 == t2); }

class AssessmentRule {
  public:
  AssessmentRule(const std::array<double,8>& g_probs) : good_probs(g_probs) {};
  // {P(G|B,B,D), P(G|B,B,C), P(G|B,G,D), P(G|B,G,C), P(G|G,B,D), P(G|G,B,C), P(G|G,G,D), P(G|G,G,C)}

  AssessmentRule(const std::map< std::tuple<Reputation,Reputation,Action>, double> & g_probs) {
    if (g_probs.size() != 8) {
      throw std::runtime_error("AssessmentRule: g_probs must have 8 elements");
    }
    good_probs[0] = g_probs.at({Reputation::B, Reputation::B, Action::D});
    good_probs[1] = g_probs.at({Reputation::B, Reputation::B, Action::C});
    good_probs[2] = g_probs.at({Reputation::B, Reputation::G, Action::D});
    good_probs[3] = g_probs.at({Reputation::B, Reputation::G, Action::C});
    good_probs[4] = g_probs.at({Reputation::G, Reputation::B, Action::D});
    good_probs[5] = g_probs.at({Reputation::G, Reputation::B, Action::C});
    good_probs[6] = g_probs.at({Reputation::G, Reputation::G, Action::D});
    good_probs[7] = g_probs.at({Reputation::G, Reputation::G, Action::C});
  }

  std::array<double,8> good_probs;

  AssessmentRule Clone() const { return AssessmentRule(good_probs); }
  AssessmentRule SwapGB() const {
    std::array<double,8> new_good_probs = {1.0-good_probs[6], 1.0-good_probs[7], 1.0-good_probs[4], 1.0-good_probs[5],
                                           1.0-good_probs[2], 1.0-good_probs[3], 1.0-good_probs[0], 1.0-good_probs[1]};
    return {new_good_probs};
  }
  double GProb(const Reputation& rep_d, const Reputation& rep_r, const Action& act) const {
    size_t idx = static_cast<size_t>(rep_d) * 4 + static_cast<size_t>(rep_r) * 2 + static_cast<size_t>(act);
    return good_probs[idx];
  }
  void SetGProb(const Reputation& rep_d, const Reputation& rep_r, const Action& act, double g_prob) {
    size_t idx = static_cast<size_t>(rep_d) * 4 + static_cast<size_t>(rep_r) * 2 + static_cast<size_t>(act);
    good_probs[idx] = g_prob;
  }

  std::string Inspect() const {
    std::stringstream ss;
    ss << "AssessmentRule: " << ID() << std::endl;
    for (size_t i = 0; i < 8; i++) {
      Reputation rep_d = static_cast<Reputation>(i / 4);
      Reputation rep_r = static_cast<Reputation>((i/2) % 2);
      Action act = static_cast<Action>(i % 2);
      ss << "(" << rep_d << "->" << rep_r << "," << act << ") : " << good_probs[i];
      if (i % 4 == 3) { ss << std::endl; }
      else { ss << "\t"; }
    }
    return ss.str();
  }

  bool IsSecondOrder() const {
    if (good_probs[0] == good_probs[4] && good_probs[1] == good_probs[5] &&
        good_probs[2] == good_probs[6] && good_probs[3] == good_probs[7]) { return true; }
    else { return false; }
  }

  bool IsDeterministic() const {
    if (( good_probs[0] == 0.0 || good_probs[0] == 1.0 ) &&
        ( good_probs[1] == 0.0 || good_probs[1] == 1.0 ) &&
        ( good_probs[2] == 0.0 || good_probs[2] == 1.0 ) &&
        ( good_probs[3] == 0.0 || good_probs[3] == 1.0 ) &&
        ( good_probs[4] == 0.0 || good_probs[4] == 1.0 ) &&
        ( good_probs[5] == 0.0 || good_probs[5] == 1.0 ) &&
        ( good_probs[6] == 0.0 || good_probs[6] == 1.0 ) &&
        ( good_probs[7] == 0.0 || good_probs[7] == 1.0 ) ) {
      return true;
    }
    return false;
  }

  int ID() const {
    if (!IsDeterministic()) { return -1; }
    int id = 0;
    for (size_t i = 0; i < 8; i++) {
      if (good_probs[i] == 1.0) { id += 1 << i; }
    }
    return id;
  }

  static AssessmentRule MakeDeterministicRule(int id) {
    if (id < 0 || id > 255) {
      throw std::runtime_error("AssessmentRuleDet: id must be between 0 and 255");
    }
    std::array<double,8> good_probs = {0.0};
    for (size_t i = 0; i < 8; i++) {
      if ((id >> i) % 2) { good_probs[i] = 1.0; }
      else { good_probs[i] = 0.0; }
    }
    return AssessmentRule{ good_probs };
  }

  static AssessmentRule AllGood() {
    return AssessmentRule::MakeDeterministicRule(0b11111111);
  }
  static AssessmentRule AllBad() {
    return AssessmentRule::MakeDeterministicRule(0b00000000);
  }
  static AssessmentRule ImageScoring() {
    return AssessmentRule{{0,1,0,1,0,1,0,1}};
  }
  static AssessmentRule KeepRecipient() {
    return AssessmentRule{{0,0,1,1,0,0,1,1}};
  }
};

bool operator==(const AssessmentRule& t1, const AssessmentRule& t2) {
  return t1.good_probs[0] == t2.good_probs[0] &&
         t1.good_probs[1] == t2.good_probs[1] &&
         t1.good_probs[2] == t2.good_probs[2] &&
         t1.good_probs[3] == t2.good_probs[3] &&
         t1.good_probs[4] == t2.good_probs[4] &&
         t1.good_probs[5] == t2.good_probs[5] &&
         t1.good_probs[6] == t2.good_probs[6] &&
         t1.good_probs[7] == t2.good_probs[7];
}
bool operator!=(const AssessmentRule& t1, const AssessmentRule& t2) {
  return !(t1 == t2);
}

// Norm3rd is a set of AssessmentRule & ActionRule
class Norm3rd {
public:
  Norm3rd(const AssessmentRule &R_d, const ActionRule &act)
      : Rd(R_d.good_probs), P(act.coop_probs) {};
  Norm3rd(const Norm3rd &rhs) : Rd(rhs.Rd), P(rhs.P) {};
  static Norm3rd FromSerialized(const std::array<double,12>& serialized) {
    std::array<double,8> Rd_serialized{};
    std::copy(serialized.begin(), serialized.begin()+8, Rd_serialized.begin());
    std::array<double,4> P_serialized{};
    std::copy(serialized.begin()+8, serialized.begin()+12, P_serialized.begin());
    return Norm3rd{AssessmentRule{Rd_serialized}, ActionRule{P_serialized}};
  }
  AssessmentRule Rd;  // assessment rules to assess donor
  ActionRule P;  // action rule
  std::string Inspect() const {
    std::stringstream ss;
    if (IsDeterministic()) {
      ss << "Norm3rd: 0x" << std::setfill('0') << std::setw(5) << std::hex << ID() << " " << std::dec << ID() << " : " << GetName() << std::endl;
      ss << std::resetiosflags(std::ios_base::fmtflags(-1));
      for (int i = 3; i >= 0; i--) {
        Reputation X = static_cast<Reputation>(i / 2);
        Reputation Y = static_cast<Reputation>(i % 2);
        Action A = (CProb(X, Y) == 1.0) ? Action::C : Action::D;
        Action notA = FlipAction(A);
        Reputation donor_rep = Rd.GProb(X, Y, A) == 1.0 ? Reputation::G : Reputation::B;
        Reputation donor_rep_not = Rd.GProb(X, Y, notA) == 1.0 ? Reputation::G : Reputation::B;
        ss << "(" << X << "->" << Y << "):" << A << donor_rep << donor_rep_not << "\t";
      }
      ss << "\n";
    }
    for (int i = 3; i >= 0; i--) {
      auto X = static_cast<Reputation>(i / 2);
      auto Y = static_cast<Reputation>(i % 2);
      double c_prob = P.CProb(X, Y);
      double gd_prob_d = Rd.GProb(X, Y, Action::D);
      double gd_prob_c = Rd.GProb(X, Y, Action::C);
      ss << std::setprecision(3) << std::fixed;
      ss << "(" << X << "->" << Y << "): "
         << "P:" << c_prob << " : "
         << "R (c:" << gd_prob_c << ",d:" << gd_prob_d << ")\n";
    }
    return ss.str();
  }
  std::string InspectComparison(const Norm3rd& other) const {
    // compare this norm and `other` norm and highlight the differences
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(5) << std::hex;
    ss << "Norm3rd: 0x" << ID() << " : " << GetName() << " vs 0x" << other.ID() << " : " << other.GetName() << std::endl;
    ss << std::resetiosflags(std::ios_base::fmtflags(-1));

    for (int i = 3; i >= 0; i--) {
      auto X = static_cast<Reputation>(i / 2);
      auto Y = static_cast<Reputation>(i % 2);
      double c_prob1 = P.CProb(X, Y);
      double gd_prob_d1 = Rd.GProb(X, Y, Action::D);
      double gd_prob_c1 = Rd.GProb(X, Y, Action::C);

      double c_prob2 = other.P.CProb(X, Y);
      double gd_prob_d2 = other.Rd.GProb(X, Y, Action::D);
      double gd_prob_c2 = other.Rd.GProb(X, Y, Action::C);
      ss << std::setprecision(2) << std::fixed;
      ss << "(" << X << "->" << Y << "): ";
      ss << "P:";
      auto highlight_diff = [&ss](double a, double b) {
        if (a != b) {
          ss << "\033[1;31m";
          ss << a << " != " << b;
          ss << "\033[0m";
        } else {
          ss << a << "        ";
        }
      };
      highlight_diff(c_prob1, c_prob2);
      ss << " : ";
      ss << "R1 (c:";
      highlight_diff(gd_prob_c1, gd_prob_c2);
      ss << ",d:";
      highlight_diff(gd_prob_d1, gd_prob_d2);
      ss << ")\n";
    }
    return ss.str();
  }
  std::array<double,12> Serialize() const {
    std::array<double,12> out{};
    for (size_t i = 0; i < 8; i++) {
      out[i] = Rd.good_probs[i];
    }
    for (size_t i = 0; i < 4; i++) {
      out[i+8] = P.coop_probs[i];
    }
    return out;
  }
  double CProb(Reputation donor, Reputation recipient) const { return P.CProb(donor, recipient); }
  double GProbDonor(Reputation donor, Reputation recipient, Action act) const {
    return Rd.GProb(donor, recipient, act);
  }
  Norm3rd SwapGB() const {
    return Norm3rd(Rd.SwapGB(), P.SwapGB());
  }
  bool IsSecondOrder() const {
    return P.IsSecondOrder() && Rd.IsSecondOrder();
  }
  bool IsDeterministic() const {
    return P.IsDeterministic() && Rd.IsDeterministic();
  }
  Norm3rd RescaleWithError(double mu_e,
                        double mu_a_donor,
                        double mu_a_recip) const { // return the norm that takes into account error probabilities
    std::array<double, 8> good_probs_donor = {0.0};
    for (size_t i = 0; i < 8; i++) {
      good_probs_donor[i] = (1.0 - 2.0 * mu_a_donor) * Rd.good_probs[i] + mu_a_donor;
    }
    return Norm3rd{{good_probs_donor}, P.RescaleWithError(mu_e)};
  }
  int ID() const {
    if (!IsDeterministic()) {
      return -1;
    }
    int id = 0;
    id += Rd.ID() << 4;
    id += P.ID();
    return id;
  }
  static Norm3rd ConstructFromID(int id) {
    if (id < 0 || id >= (1 << 12)) {
      throw std::runtime_error("Norm3rd: id must be between 0 and 2^12-1");
    }
    int Rd_id = (id >> 4) & 0xFF;
    int P_id = id & 0xF;
    return Norm3rd(AssessmentRule::MakeDeterministicRule(Rd_id),
                   ActionRule::MakeDeterministicRule(P_id));
  }
  static Norm3rd AllC() {
    return Norm3rd(AssessmentRule{{1, 1, 1, 1, 1, 1, 1, 1}},
                   ActionRule{{1, 1, 1, 1}});
  }
  static Norm3rd AllD() {
    return Norm3rd(AssessmentRule{{0, 0, 0, 0, 0, 0, 0, 0}},
                   ActionRule{{0, 0, 0, 0}});
  }
  static Norm3rd AllG() {
    // Always assess G, but action rule is Discriminator
    // it may defect under assessment error
    return Norm3rd(AssessmentRule::AllGood(),
                   ActionRule::DISC());
  }
  static Norm3rd AllB() {
    // Always assess B, but action rule is Discriminator
    // it may defect under assessment error
    return Norm3rd(AssessmentRule::AllBad(),
                   ActionRule::DISC());
  }
  static Norm3rd ImageScoring() {
    return Norm3rd(AssessmentRule::ImageScoring(),
                   ActionRule::DISC());
  }
  static Norm3rd Random(double p = 0.5) {
    return Norm3rd{ AssessmentRule::AllBad(),
                    ActionRule{{p, p, p, p}} };
  }
  static Norm3rd L1() {
    return Norm3rd({{0, 1, 0, 1, 1, 1, 0, 1}},
                   {{1, 1, 0, 1}});
  }
  static Norm3rd L2() {
    return Norm3rd({{0, 1, 0, 1, 1, 0, 0, 1}},
                   {{1, 1, 0, 1}});
  }
  static Norm3rd L3() {
    return Norm3rd({{1, 1, 0, 1, 1, 1, 0, 1}},
                   {{0, 1, 0, 1}});
  }
  static Norm3rd L4() {
    return Norm3rd({{1, 0, 0, 1, 1, 1, 0, 1}},
                   {{0, 1, 0, 1}});
  }
  static Norm3rd L5() {
    return Norm3rd({{1, 1, 0, 1, 1, 0, 0, 1}},
                   {{0, 1, 0, 1}});
  }
  static Norm3rd L6() {
    return Norm3rd({{1, 0, 0, 1, 1, 0, 0, 1}},
                   {{0, 1, 0, 1}});
  }
  static Norm3rd L7() {
    return Norm3rd({{0, 0, 0, 1, 1, 1, 0, 1}},
                   {{0, 1, 0, 1}});
  }
  static Norm3rd L8() {
    return Norm3rd({{0, 0, 0, 1, 1, 0, 0, 1}},
                   {{0, 1, 0, 1}});
  }
  static Norm3rd SecondarySixteen(int i) {
    if (i <= 0 || i > 16) { throw std::runtime_error("Norm3rd: i must be between 1 and 16"); }
    double R_GB_C = static_cast<double>( ((i - 1) >> 3) & 0b1 );
    double R_BG_C = static_cast<double>( ((i - 1) >> 2) & 0b1 );
    double R_BB_C = static_cast<double>( ((i - 1) >> 1) & 0b1 );
    double R_BB_D = static_cast<double>( ((i - 1) >> 0) & 0b1 );
    double P_BB = 0.0;
    if (R_BB_C == 1.0 && R_BB_D == 0.0) { P_BB = 1.0; }
    return Norm3rd({{R_BB_D, R_BB_C, 1, R_BG_C, 1, R_GB_C, 0, 1}},
                   {{P_BB, 0, 0, 1}});
  }
  static Norm3rd GenerousScoring(double benefit = 5.0) {
    double q = 1.0 / benefit;
    AssessmentRule gsco(std::array<double,8>{{q, 1.0, q, 1.0, q, 1.0, q, 1.0}});
    return Norm3rd(gsco, ActionRule::DISC());
  }
  bool IsGenerousScoring() const {
    if( P != ActionRule::DISC() ) {
      return false;
    }
    if (Rd.good_probs[1] < 1.0 || Rd.good_probs[3] < 1.0 || Rd.good_probs[5] < 1.0 || Rd.good_probs[7] < 1.0) {
      return false;
    }
    if (Rd.good_probs[0] != Rd.good_probs[2] || Rd.good_probs[0] != Rd.good_probs[4] || Rd.good_probs[0] != Rd.good_probs[6]) {
      return false;
    }
    return true;
  }

  static const std::vector<std::pair<int, std::string> > NormNames;
  std::string GetName() const {
    if (IsDeterministic()) {
      for (auto &p : NormNames) {
        if (p.first == ID() || p.first == SwapGB().ID() ) {
          return p.second;
        }
      }
      if (P.ID() == 0) { return "AllD"; }
      if (P.ID() == 15) { return "AllC"; }
    }
    if (IsGenerousScoring()) {
      double benefit = 1.0 / Rd.good_probs[0];
      std::ostringstream oss;
      oss << "GSCO-" << std::fixed << std::setprecision(1) << benefit;
      return oss.str();
    }
    return "";
  }
  static Norm3rd ConstructFromName(const std::string &name) {
    auto found = std::find_if(Norm3rd::NormNames.begin(), Norm3rd::NormNames.end(), [&](auto &s) {
      return s.second == name;
    });

    std::regex re_gsco(R"(^GSCO-(\d+(?:\.\d+)?)$)");
    std::regex re_random(R"(^RANDOM-(0(\.\d+)?|1(\.0+)?)$)");  // RANDOM-0.5
    std::smatch matches;

    if (found != Norm3rd::NormNames.end()) {
      return Norm3rd::ConstructFromID(found->first);
    }
    else if (std::regex_match(name, matches, re_gsco)) {
      double benefit = std::stod(matches[1]);
      return Norm3rd::GenerousScoring(benefit);
    }
    else if (std::regex_match(name, matches, re_random)) {
      double p = std::stod(matches[1]);
      return Norm3rd::Random(p);
    }
    else {
      std::cerr << "Unknown norm: " << name << std::endl;
      std::cerr << "Available norms are: " << std::endl;
      for (auto &n : Norm3rd::NormNames) {
        std::cerr << "\t" << n.second << std::endl;
      }
      std::cerr << "\t" << "GSCO-benefit" << std::endl;
      throw std::runtime_error("Unknown norm");
    }
  }

  static Norm3rd ParseNormString(const std::string& str, bool swap_gb = false) {
    std::regex re_d(R"(\d+)"); // regex for digits
    std::regex re_x(R"(^0x[0-9a-fA-F]+$)");  // regex for digits in hexadecimal
    // regular expression for 12 floating point numbers separated by space
    std::regex re_a(R"(^(0|1)(\.\d+)?( (0|1)(\.\d+)?){11}$)");
    Norm3rd norm = Norm3rd::AllC();
    if (std::regex_match(str, re_d)) {
      int id = std::stoi(str);
      norm = Norm3rd::ConstructFromID(id);
    }
    else if (std::regex_match(str, re_x)) {
      int id = std::stoi(str, nullptr, 16);
      norm = Norm3rd::ConstructFromID(id);
    }
    else if (std::regex_match(str, re_a)) {
      std::istringstream iss(str);
      std::array<double,12> serialized = {};
      for (int i = 0; i < 12; ++i) {
        iss >> serialized[i];
      }
      norm = Norm3rd::FromSerialized(serialized);
    }
    else {
      norm = Norm3rd::ConstructFromName(str);
    }
    if (swap_gb) {
      norm = norm.SwapGB();
    }
    return norm;
  }
};

const std::vector<std::pair<int,std::string> > Norm3rd::NormNames = {{
                                                                    {AllC().ID(), "AllC"},
                                                                    {AllD().ID(), "AllD"},
                                                                    {AllG().ID(), "AllG"},
                                                                    {AllB().ID(), "AllB"},
                                                                    {ImageScoring().ID(), "ImageScoring"},
                                                                    {L1().ID(), "L1"},
                                                                    {L2().ID(), "L2"},
                                                                    {L3().ID(), "L3"},
                                                                    {L4().ID(), "L4"},
                                                                    {L5().ID(), "L5"},
                                                                    {L6().ID(), "L6"},
                                                                    {L7().ID(), "L7"},
                                                                    {L8().ID(), "L8"},
                                                                    {SecondarySixteen(1).ID(), "S1"},
                                                                    {SecondarySixteen(2).ID(), "S2"},
                                                                    {SecondarySixteen(3).ID(), "S3"},
                                                                    {SecondarySixteen(4).ID(), "S4"},
                                                                    {SecondarySixteen(5).ID(), "S5"},
                                                                    {SecondarySixteen(6).ID(), "S6"},
                                                                    {SecondarySixteen(7).ID(), "S7"},
                                                                    {SecondarySixteen(8).ID(), "S8"},
                                                                    {SecondarySixteen(9).ID(), "S9"},
                                                                    {SecondarySixteen(10).ID(), "S10"},
                                                                    {SecondarySixteen(11).ID(), "S11"},
                                                                    {SecondarySixteen(12).ID(), "S12"},
                                                                    {SecondarySixteen(13).ID(), "S13"},
                                                                    {SecondarySixteen(14).ID(), "S14"},
                                                                    {SecondarySixteen(15).ID(), "S15"},
                                                                    {SecondarySixteen(16).ID(), "S16"},
                                                                }};
bool operator==(const Norm3rd& n1, const Norm3rd& n2) {
  return n1.P == n2.P && n1.Rd == n2.Rd;
}
bool operator!=(const Norm3rd& n1, const Norm3rd& n2) {
  return !(n1 == n2);
}

#endif
