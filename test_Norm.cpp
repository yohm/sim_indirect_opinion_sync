#include <iostream>
#include <bitset>
#include <regex>
#include <vector>
#include <set>
#include "Norm.hpp"

#include <gtest/gtest.h>

constexpr Reputation B = Reputation::B, G = Reputation::G;
constexpr Action C = Action::C, D = Action::D;

TEST(ActionRule, action_rules) {
  ActionRule P(0.9, 0.1);
  std::cerr << P.Inspect();

  EXPECT_DOUBLE_EQ( P.CProb(G), 0.9 );
  EXPECT_DOUBLE_EQ( P.CProb(B), 0.1 );

  P.SetCProb(G, 0.3);
  EXPECT_DOUBLE_EQ( P.CProb(G), 0.3 );

  ActionRule P2 = P.SwapGB();
  EXPECT_DOUBLE_EQ( P2.CProb(G), 0.1 );
  EXPECT_DOUBLE_EQ( P2.CProb(B), 0.3 );

  ActionRule P3(1.0, 0.0);
  std::cerr << P3.Inspect();
  EXPECT_TRUE( P3.IsDeterministic() );
  EXPECT_EQ( P3.ID(), 0b10 );

  auto disc = ActionRule::DISC();
  EXPECT_TRUE( disc.IsDeterministic() );
  EXPECT_EQ( disc.ID(), 0b10 );

  auto adisc = ActionRule::ADISC();
  EXPECT_TRUE( adisc.IsDeterministic() );
  EXPECT_EQ( adisc.ID(), 0b01 );

  auto allc = ActionRule::ALLC();
  EXPECT_TRUE( allc.IsDeterministic() );
  EXPECT_EQ( allc.ID(), 0b11 );

  auto alld = ActionRule::ALLD();
  EXPECT_TRUE( alld.IsDeterministic() );
  EXPECT_EQ( alld.ID(), 0b00 );

  auto P4 = ActionRule::MakeDeterministicRule(2);
  EXPECT_TRUE( P4.IsDeterministic() );
  EXPECT_EQ( P4.ID(), 2 );
}

TEST(AssessmentRule, assessment_rules) {
  AssessmentRule R(0.8, 0.6, 0.4, 0.2);
  std::cout << R.Inspect();

  EXPECT_DOUBLE_EQ( R.GProb(B, D), 0.2 );
  EXPECT_DOUBLE_EQ( R.GProb(B, C), 0.4 );
  EXPECT_DOUBLE_EQ( R.GProb(G, D), 0.6 );
  EXPECT_DOUBLE_EQ( R.GProb(G, C), 0.8 );

  R.SetGProb(G, C, 0.9);
  EXPECT_DOUBLE_EQ( R.GProb(G, C), 0.9 );

  AssessmentRule R2 = R.SwapGB();
  EXPECT_DOUBLE_EQ( R2.GProb(G, C), 1.0 - R.GProb(B, C) );

  AssessmentRule Q3(1.0, 0.0, 1.0, 0.0);
  EXPECT_TRUE( Q3.IsDeterministic() );
  EXPECT_EQ( Q3.ID(), 0b1010 );

  auto allg = AssessmentRule::AllGood();
  EXPECT_TRUE( allg.IsDeterministic() );
  EXPECT_EQ( allg.ID(), 0b1111 );
  auto allb = AssessmentRule::AllBad();
  EXPECT_TRUE( allb.IsDeterministic() );
  EXPECT_EQ( allb.ID(), 0b0000 );
  auto is = AssessmentRule::ImageScoring();
  EXPECT_TRUE( is.IsDeterministic() );
  EXPECT_EQ( is.ID(), 0b1010 );
}

TEST(Norm, base) {
  Norm n{
      AssessmentRule{0.8, 0.6, 0.4, 0.2},
      ActionRule{0.9, 0.1}
  };
  EXPECT_FALSE(n.IsDeterministic());
  EXPECT_DOUBLE_EQ(n.Gprob(G, C), 0.8);
  EXPECT_DOUBLE_EQ(n.CProb(G), 0.9);
  // std::cout << n.Inspect();
}

TEST(Norm, ImageScoring) {
  Norm is = Norm::ImageScoring();
  EXPECT_TRUE( is.IsDeterministic() );
  EXPECT_EQ( is.ID(), 0b1010'10);
}

TEST(Norm, AllG) {
  Norm allg = Norm::AllG();
  EXPECT_TRUE( allg.IsDeterministic() );
  EXPECT_EQ( allg.ID(), 0b1111'10);
}

TEST(Norm, AllB) {
  Norm allb = Norm::AllB();
  EXPECT_TRUE( allb.IsDeterministic() );
  EXPECT_EQ( allb.ID(), 0b0000'10);
}

TEST(Norm, LeadingEight) {
  std::array<Norm, 2> leading_eight = {Norm::L3(), Norm::L6()};
  std::array<int, 2> l8_ids = {0b1011'10, 0b1001'10 };
  std::array<std::string,2> names = {"L3", "L6"};
  for (size_t i = 0; i < 2; i++) {
    auto l = leading_eight[i];
    EXPECT_TRUE( l.IsDeterministic() );

    EXPECT_EQ( l.GetName(), names[i]);
    EXPECT_EQ( l.ID(), l8_ids[i] );
    EXPECT_EQ( Norm::ConstructFromID(l.ID()), l );

    EXPECT_DOUBLE_EQ( l.CProb(G), 1.0 );
    EXPECT_DOUBLE_EQ( l.CProb(B), 0.0 );
    EXPECT_DOUBLE_EQ(l.Gprob(G, C), 1.0 );
    EXPECT_DOUBLE_EQ(l.Gprob(G, D), 0.0 );  // identification of defectors
    EXPECT_DOUBLE_EQ(l.Gprob(B, D), 1.0 );  // justified punishment
  }
}

TEST(Norm, rescaling) {
  Norm n = Norm::L6();
  auto rescaled = n.RescaleWithError(0.1, 0.02);
  EXPECT_DOUBLE_EQ(rescaled.CProb(G), 0.9);
  EXPECT_DOUBLE_EQ(rescaled.CProb(B), 0.0);

  EXPECT_DOUBLE_EQ(rescaled.Gprob(G, C), 0.98);
  EXPECT_DOUBLE_EQ(rescaled.Gprob(G, D), 0.02);
  EXPECT_DOUBLE_EQ(rescaled.Gprob(B, C), 0.02);
  EXPECT_DOUBLE_EQ(rescaled.Gprob(B, D), 0.98);
}

TEST(Norm, GenerousScoring) {
  Norm gsco = Norm::GenerousScoring(5.0);
  EXPECT_TRUE( gsco.IsGenerousScoring() );
  std::string name = gsco.GetName();
  EXPECT_EQ( name, "GSCO-5.0" );

  gsco = Norm::GenerousScoring(1.5);
  EXPECT_EQ( gsco.GetName(), "GSCO-1.5" );

  Norm n = Norm::ConstructFromName("GSCO-1.5");
  EXPECT_TRUE( n.IsGenerousScoring() );
  EXPECT_EQ(n.GetName(), "GSCO-1.5");
  EXPECT_DOUBLE_EQ(n.Gprob(G, D), 2.0 / 3.0 );
}

TEST(Norm, Deterministic2ndOrder) {
  std::vector<Norm> norms = Norm::DeterministicNorms();
  std::set<int> ids;
  for (const Norm& n: norms) {
    EXPECT_TRUE(n.IsDeterministic());
    ids.insert(n.ID());
  }
  EXPECT_EQ( norms.size(), 36 );
  EXPECT_EQ( ids.size(), 36 );
}

TEST(Norm, ParseNormString) {
  EXPECT_EQ( Norm::ParseNormString("AllG").GetName(), "AllG" );
  EXPECT_EQ( Norm::ParseNormString("L3").GetName(), "L3" );

  EXPECT_EQ( Norm::ParseNormString("GSCO-1.5").GetName(), "GSCO-1.5" );
  EXPECT_EQ( Norm::ParseNormString("33").ID(), 33 );

  EXPECT_EQ( Norm::ParseNormString("SECOND-0").ID(), 0b0000'10 );
  EXPECT_EQ( Norm::ParseNormString("SECOND-3").ID(), 0b0011'10 );
  EXPECT_EQ( Norm::ParseNormString("SECOND-9").ID(), 0b1001'10 );
  EXPECT_EQ( Norm::ParseNormString("SECOND-9").GetName(), "L6" );
  EXPECT_EQ( Norm::ParseNormString("SECOND-11").ID(), 0b1011'10 );
  EXPECT_EQ( Norm::ParseNormString("SECOND-11").GetName(), "L3" );
  EXPECT_EQ( Norm::ParseNormString("SECOND-15").ID(), 0b1111'10 );

  {
    double p = 0.7;
    auto n = Norm::ParseNormString("RANDOM-0.7");
    EXPECT_EQ( n.CProb(G), p );
    EXPECT_EQ( n.CProb(B), p );
  }
  {
    double p = 1.0;
    auto n = Norm::ParseNormString("RANDOM-1.0");
    EXPECT_EQ( n.CProb(G), p );
    EXPECT_EQ( n.CProb(B), p );
  }

  const std::string s = "0.8 0.6 0.4 0.2 1.0 0.0";
  Norm n = Norm::ParseNormString(s);
  EXPECT_DOUBLE_EQ( n.CProb(G), 1.0 );
  EXPECT_DOUBLE_EQ( n.CProb(B), 0.0 );
  EXPECT_DOUBLE_EQ(n.Gprob(G, C), 0.8 );
  EXPECT_DOUBLE_EQ(n.Gprob(G, D), 0.6 );
  EXPECT_DOUBLE_EQ(n.Gprob(B, C), 0.4 );
  EXPECT_DOUBLE_EQ(n.Gprob(B, D), 0.2 );
}
