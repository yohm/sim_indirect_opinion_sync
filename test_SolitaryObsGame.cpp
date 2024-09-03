#include <iostream>
#include <fstream>
#include <chrono>
#include <regex>
#include "SolitaryObsGame.hpp"

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
  double h = SolitaryObsGame::HscoreMonomorphic(norm);
  EXPECT_NEAR(h, expected_good_rep, 0.01);
  double pc = SolitaryObsGame::SelfCooperationLevel(norm, h);
  EXPECT_NEAR(pc, expected_c_level, 0.01);
}

TEST_F(TestFixture, BasicNorms) {
  auto random = Norm{0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
  test_SelfCooperationLevel(random, 0.5, 0.5);
  {
    // irrespective of initial h, it converges to 0.5
    std::vector<double> init_h = {0.5, 0.1, 0.9};
    for (double h: init_h) {
      auto ht = SolitaryObsGame::CalculateHDynamics(random, h, 0, 0.01, 1000);
      std::cerr << ht[0] << std::endl;
      EXPECT_NEAR(ht[0], 0.5, 0.01);
    }
  }

  test_SelfCooperationLevel(Norm::AllC(), 1.0, 1.0);
  test_SelfCooperationLevel(Norm::AllD(), 0.0, 0.0);

  test_SelfCooperationLevel(Norm::ImageScoring().RescaleWithError(0.0, 0.01), 0.5, 0.5);
}

TEST_F(TestFixture, PolymorphicPopulations) {
  auto random = Norm{0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
  auto alld = Norm::AllD();
  auto h_timeseries = SolitaryObsGame::CalculateHDynamicsPolymorphic(random, alld, 0.5);
  EXPECT_EQ(h_timeseries.size(), 1);
  auto [h_rr, h_rm, h_mr, h_mm] = h_timeseries.at(0);
  EXPECT_NEAR(h_rr, 0.5, 0.01);
  EXPECT_NEAR(h_rm, 0.5, 0.01);
  EXPECT_NEAR(h_mr, 0.0, 0.01);
  EXPECT_NEAR(h_mm, 0.0, 0.01);
  SolitaryObsGame::CooperationLevelsPolymorphic(random, alld, 0.5, h_timeseries[0]);

  auto [pc_rr,pc_rm,pc_mr,pc_mm] = SolitaryObsGame::CooperationLevelsPolymorphic(random, alld, 0.5, h_timeseries[0]);
  EXPECT_NEAR(pc_rr, 0.5, 0.01);
  EXPECT_NEAR(pc_rm, 0.5, 0.01);
  EXPECT_NEAR(pc_mr, 0.0, 0.01);
  EXPECT_NEAR(pc_mm, 0.0, 0.01);

  auto allc = Norm::AllC();
  auto is = Norm::ImageScoring().RescaleWithError(0.0, 0.01);
  h_timeseries = SolitaryObsGame::CalculateHDynamicsPolymorphic(allc, is, 0.5);
  EXPECT_EQ(h_timeseries.size(), 1);
  std::tie(h_rr, h_rm, h_mr, h_mm) = h_timeseries.at(0);
  EXPECT_NEAR(h_rr, 1.0, 0.01);
  EXPECT_NEAR(h_rm, 1.0, 0.01);
  EXPECT_NEAR(h_mr, 1.0, 0.02);
  EXPECT_NEAR(h_mm, 1.0, 0.04);
}
