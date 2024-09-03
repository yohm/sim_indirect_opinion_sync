#include <iostream>
#include <fstream>
#include <chrono>
#include <regex>
#include "PrivRepGame.hpp"

#include <gtest/gtest.h>

class TestFixture : public ::testing::Test {
  protected:
  void SetUp() override {
    start_time_ = std::chrono::high_resolution_clock::now();
  }

  void TearDown() override {
    auto end_time = std::chrono::high_resolution_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time_);
    std::cout << "Test " << ::testing::UnitTest::GetInstance()->current_test_info()->name() << " elapsed time: " << elapsed_time.count() << " ms" << std::endl;
  }

  private:
  std::chrono::time_point<std::chrono::high_resolution_clock> start_time_;
};


void test_SelfCooperationLevel(const Norm& norm, double expected_c_level, double expected_good_rep) {
  PrivateRepGame priv_game( {{norm, 50}}, 123456789ull);
  priv_game.Update(1e3, 0.9, 0.0, 0.02, 0.0, false);
  priv_game.ResetCounts();
  priv_game.Update(1e3, 0.9, 0.0, 0.02, 0.0, true);
  // IC( priv_game.NormCooperationLevels(), priv_game.NormAverageReputation() );
  EXPECT_NEAR( priv_game.SystemWideCooperationLevel(), expected_c_level, 0.02);
  EXPECT_NEAR( priv_game.NormCooperationLevels()[0][0], expected_c_level, 0.02);
  EXPECT_NEAR( priv_game.NormAverageReputation()[0][0], expected_good_rep, 0.02);
}

TEST_F(TestFixture, RandomNorm) {
  test_SelfCooperationLevel(Norm::Random(), 0.5, 0.02);
}

TEST_F(TestFixture, LeadingEight) {
  test_SelfCooperationLevel(Norm::L3(), 0.95, 0.95);
  test_SelfCooperationLevel(Norm::L6(), 0.50, 0.50);
}

TEST_F(TestFixture, RandomNonIdenticalPermutations) {
  for (int t = 0; t < 100; t++) {
    std::mt19937_64 rng(123456789ull + t);
    auto [perm1, perm2] = PrivateRepGame::RandomNonIdenticalPermutations(2, rng);
    EXPECT_EQ(perm1.size(), 2);
    EXPECT_EQ(perm2.size(), 2);
    std::set<size_t> s1, s2;
    for (int i = 0; i < perm1.size(); i++) {
      EXPECT_NE(perm1[i], perm2[i]);
      s1.insert(perm1[i]);
      s2.insert(perm2[i]);
    }
    EXPECT_EQ(s1.size(), 2);
    EXPECT_EQ(s2.size(), 2);
  }
  for (int t = 0; t < 100; t++) {
    std::mt19937_64 rng(123456789ull + t);
    auto [perm1, perm2] = PrivateRepGame::RandomNonIdenticalPermutations(3, rng);
    EXPECT_EQ(perm1.size(), 3);
    EXPECT_EQ(perm2.size(), 3);
    std::set<size_t> s1, s2;
    for (int i = 0; i < perm1.size(); i++) {
      EXPECT_NE(perm1[i], perm2[i]);
      s1.insert(perm1[i]);
      s2.insert(perm2[i]);
    }
    EXPECT_EQ(s1.size(), 3);
    EXPECT_EQ(s2.size(), 3);
  }
  for (int t = 0; t < 100; t++) {
    std::mt19937_64 rng(123456789ull + t);
    auto [perm1, perm2] = PrivateRepGame::RandomNonIdenticalPermutations(50, rng);
    EXPECT_EQ(perm1.size(), 50);
    EXPECT_EQ(perm2.size(), 50);
    std::set<size_t> s1, s2;
    for (int i = 0; i < perm1.size(); i++) {
      EXPECT_NE(perm1[i], perm2[i]);
      s1.insert(perm1[i]);
      s2.insert(perm2[i]);
    }
    EXPECT_EQ(s1.size(), 50);
    EXPECT_EQ(s2.size(), 50);
  }
}

TEST(Reputation, h_scores) {
  {
    Norm norm = Norm::AllD();
    PrivateRepGame priv_game({{norm, 50}}, 123456789ull);
    priv_game.Update(1e3, 1.0, 0.0, 0.02, 0.0, false);
    priv_game.ResetCounts();
    priv_game.Update(1e3, 1.0, 0.0, 0.02, 0.0, true);

    auto gi_vec = priv_game.NormAverageGoodness();
    double h = gi_vec[0].first;
    double hg = (50.0 - 1.0) / (50.0 - 2.0) * gi_vec[0].second / h - 1.0 / (50 - 2);
    std::cerr << h << " " << hg << std::endl;
    EXPECT_NEAR(h, 0.02, 0.01);
    EXPECT_NEAR(hg, 0.02, 0.01);
  }
  {
    Norm norm = Norm::AllC();
    PrivateRepGame priv_game({{norm, 50}}, 123456789ull);
    priv_game.Update(1e3, 1.0, 0.0, 0.02, 0.0, false);
    priv_game.ResetCounts();
    priv_game.Update(1e3, 1.0, 0.0, 0.02, 0.0, true);

    auto gi_vec = priv_game.NormAverageGoodness();
    double h = gi_vec[0].first;
    double hg = (50.0 - 1.0) / (50.0 - 2.0) * gi_vec[0].second / h - 1.0 / (50 - 2);
    std::cerr << h << " " << hg << std::endl;
    EXPECT_NEAR(h, 0.98, 0.01);
    EXPECT_NEAR(hg, 0.98, 0.01);
  }
  {
    Norm norm = Norm::L3();
    PrivateRepGame priv_game({{norm, 100}}, 123456789ull);
    priv_game.Update(1e3, 1.0, 0.0, 0.02, 0.0, false);
    priv_game.ResetCounts();
    priv_game.Update(1e3, 1.0, 0.0, 0.02, 0.0, true);

    auto rep = priv_game.NormAverageReputation();
    std::cerr << rep[0][0] << std::endl;
    auto gi_vec = priv_game.NormAverageGoodness();
    double h = gi_vec[0].first;
    double hg = (50.0 - 1.0) / (50.0 - 2.0) * gi_vec[0].second / h - 1.0 / (50 - 2);
    std::cerr << h << " " << hg << std::endl;
    EXPECT_NEAR(h, 0.96, 0.01);
    EXPECT_NEAR(hg, 0.978, 0.01);
    EXPECT_GT(hg, h);

    Norm r_l3 = norm.RescaleWithError(0.0, 0.02);
    constexpr Reputation G = Reputation::G, B = Reputation::B;
    auto PR = [r_l3](Reputation m_AB, Reputation m_CB) -> double {
      constexpr Action C = Action::C, D = Action::D;
      double pc = r_l3.CProb(m_AB);
      std::cerr << "pc: " << pc << std::endl;
      return pc * r_l3.Gprob(m_CB, C) + (1.0 - pc) * r_l3.Gprob(m_CB, D);
    };
    double pr_GG = PR(G, G), pr_GB = PR(G, B), pr_BG = PR(B, G), pr_BB = PR(B, B);
    std::cerr << "pr_GG: " << pr_GG << " pr_GB: " << pr_GB << " pr_BG: " << pr_BG << " pr_BB: " << pr_BB << std::endl;
    double h_exp = h * hg * pr_GG + h * (1.0 - hg) * (pr_GB + pr_BG) + (1.0 - 2.0 * h + h * hg) * pr_BB;
    std::cerr << "h: " << h << " hg: " << hg << " h_exp: " << h_exp << std::endl;
    EXPECT_NEAR(h_exp, h, 0.01);

    double bc_upper = (1.0 - h) / (h * (1.0 - hg) * (1.0 - 0.04));
    std::cerr << "bc_upper: " << bc_upper << std::endl;
    EXPECT_NEAR(bc_upper, 2.0, 0.1);  // According to Fujimoto & Ohtsuki, b/c < 2 is the stability condition against ALLC
  }
  {
    Norm norm = Norm::L3();
    PrivateRepGame priv_game({{norm, 50}}, 123456789ull);
    priv_game.Update(1e3, 0.0, 0.0, 0.02, 0.0, false);
    priv_game.ResetCounts();
    priv_game.Update(1e3, 0.0, 0.0, 0.02, 0.0, true);

    auto rep = priv_game.NormAverageReputation();
    std::cerr << rep[0][0] << std::endl;
    auto gi_vec = priv_game.NormAverageGoodness();
    double h = gi_vec[0].first;
    double hg = (50.0 - 1.0) / (50.0 - 2.0) * gi_vec[0].second / h - 1.0 / (50 - 2);
    std::cerr << h << " " << hg << std::endl;
    EXPECT_NEAR(h, 0.87, 0.01);
    EXPECT_NEAR(hg, h, 0.01);
  }

}

TEST(Gossiping, h_scores) {
  {
    Norm norm = Norm::L6();
    double tau = 0.4;
    double mu_a = 0.02;
    PrivateRepGame priv_game({{norm, 50}}, 123456789ull);
    priv_game.Update(1e3, 0.0, 0.00, mu_a, tau, false);
    priv_game.ResetCounts();
    priv_game.Update(1e3, 0.0, 0.00, mu_a, tau, true);

    auto gi_vec = priv_game.NormAverageGoodness();
    double h = gi_vec[0].first;
    double hg = (50.0 - 1.0) / (50.0 - 2.0) * gi_vec[0].second / h - 1.0 / (50 - 2);
    double pc = priv_game.NormCooperationLevels()[0][0];
    std::cerr << h << " " << hg << " " << pc << std::endl;
    double RP_GG = 1.0 - mu_a, RP_GB = mu_a, RP_BG = mu_a, RP_BB = 1.0 - mu_a;
    double h_exp = h * hg * RP_GG + h*(1-hg) * (RP_GB + RP_BG) + (1.0-2.0*h+h*hg) * RP_BB;
    EXPECT_NEAR(h, h_exp, 0.01);
  }
  {
    Norm norm = Norm::L6();
    double tau = 0.8;
    double mu_a = 0.02;
    PrivateRepGame priv_game({{norm, 50}}, 123456789ull);
    priv_game.Update(1e3, 0.0, 0.0, mu_a, tau, false);
    priv_game.ResetCounts();
    priv_game.Update(1e3, 0.0, 0.0, mu_a, tau, true);

    auto gi_vec = priv_game.NormAverageGoodness();
    double h = gi_vec[0].first;
    double hg = (50.0 - 1.0) / (50.0 - 2.0) * gi_vec[0].second / h - 1.0 / (50 - 2);
    double pc = priv_game.NormCooperationLevels()[0][0];
    std::cerr << h << " " << hg << " " << pc << std::endl;
    // [WIP] implement the test case with the expected values
    double RP_GG = 1.0 - mu_a, RP_GB = mu_a, RP_BG = mu_a, RP_BB = 1.0 - mu_a;
    double h_exp = h * hg * RP_GG + h*(1-hg) * (RP_GB + RP_BG) + (1.0-2.0*h+h*hg) * RP_BB;
    EXPECT_NEAR(h, h_exp, 0.01);
  }
}
