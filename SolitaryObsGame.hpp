#ifndef SOLITARY_OBS_GAME_HPP
#define SOLITARY_OBS_GAME_HPP

#include <iostream>
#include <vector>
#include <array>
#include <optional>
#include <cassert>
#include "Norm.hpp"


class SolitaryObsGame {
  // mean-field analysis of the indirect reciprocity game with solitary observer
  // we do not conduct simulation, but calculate the expected cooperation level numerically
  // when you introduce the errors, please rescale the norms with the error rates
public:
  static double HscoreMonomorphic(const Norm& norm) {
    // \dot{h} = h^2 RP(GG) + h(1-h) {RP(GB) + RP(BG)} + (1-h)^2 RP(BB) - h
    //         = a h^2 + b h + c

    auto [RP_GG, RP_GB, RP_BG, RP_BB] = RP(norm.P, norm.R);

    double a = RP_GG - RP_GB - RP_BG + RP_BB;
    double b = RP_GB + RP_BG - 2.0 * RP_BB - 1.0;
    double c = RP_BB;
    // std::cerr << a << ' ' << b << ' ' << c << std::endl;
    // h^{\ast} is the solution of the quadratic equation a h^2 + b h + c = 0
    // h^{\ast} = (-b - \sqrt{b^2 - 4ac}) / 2a (when a != 0)
    // h^{\ast} = -c / b (when a == 0)
    // Although analytic expression is available, it is numerically unstable when a is small.
    // Use newton's method instead.
    const double tolerance = 2.0e-3;
    double h = (std::abs(a) > tolerance) ? ((-b - std::sqrt(b*b - 4.0*a*c)) / (2.0*a)) : -c / b;
    // solve f(h) = ah^2 + bh + c = 0 by Newton's method
    if (h > 1.0) h = 1.0;
    if (h < 0.0) h = 0.0;
    double f = a*h*h + b*h + c;
    double df = 2.0*a*h + b;
    while (std::abs(f) > 1.0e-12) {
      h -= f / df;
      f = a*h*h + b*h + c;
      df = 2.0*a*h + b;
    }
    return h;
  }

  static std::array<double,4> RP(const ActionRule& P, const AssessmentRule& R) {
    constexpr Reputation B = Reputation::B, G = Reputation::G;
    constexpr Action C = Action::C, D = Action::D;
    double PG = P.CProb(G), PB = P.CProb(B);
    double RGC = R.GProb(G, C), RGD = R.GProb(G, D);
    double RBC = R.GProb(B, C), RBD = R.GProb(B, D);
    double PR_GG = PG * RGC + (1.0 - PG) * RGD;
    double PR_GB = PG * RBC + (1.0 - PG) * RBD;
    double PR_BG = PB * RGC + (1.0 - PB) * RGD;
    double PR_BB = PB * RBC + (1.0 - PB) * RBD;
    return {PR_GG, PR_GB, PR_BG, PR_BB};
  }

  static std::vector<double> CalculateHDynamics(const Norm& norm, double h_init = 0.5, size_t output_interval = 0, double dt = 0.01, size_t n_iter = 1000) {
    auto [RP_GG, RP_GB, RP_BG, RP_BB] = RP(norm.P, norm.R);
    // Solve time evolution by Runge-Kutta method
    auto hdot = [RP_GG,RP_GB,RP_BG,RP_BB](const std::vector<double>& hv, std::vector<double>& dhv) {
      double h = hv[0];
      dhv[0] = h * h * RP_GG + h * (1.0 - h) * RP_GB + h * (1.0 - h) * RP_BG + (1.0 - h) * (1.0 - h) * RP_BB - h;
    };
    std::vector<double> hvec = {h_init};
    // \dot{h} = h^2 * PR(GG) + h(1-h) * PR(GB) + h(1-h) * PR(BG) + (1-h)^2 * RP(BB)
    // RP(XY) = P(X) * R(Y, C) + (1-P(X)) * R(Y, D)
    if (output_interval == 0) {
      hvec = SolveByRungeKutta(hdot, hvec, dt, n_iter);
      return {hvec[0]};
    }
    else {
      std::vector<double> ans;
      for (size_t t = 0; t < n_iter; t+=output_interval) {
        hvec = SolveByRungeKutta(hdot, hvec, dt, output_interval);
        ans.push_back(hvec[0]);
      }
      return ans;
    }
  }

  static double SelfCooperationLevel(const Norm& norm, std::optional<double> h = std::nullopt) {
    if (!h.has_value()) {
      h = HscoreMonomorphic(norm);
    }
    double pc = h.value() * norm.CProb(Reputation::G) + (1.0 - h.value()) * norm.CProb(Reputation::B);
    return pc;
  }

  static std::array<double,4> HscorePolymorphic(const Norm& res, const Norm& mut, double nu) {
    // res: resident norm
    // mut: mutant norm
    // nu: fraction of mutants (0, 1)
    // return: [h_rr, h_rm, h_mr, h_mm]
    auto ht = CalculateHDynamicsPolymorphic(res, mut, nu);
    return ht[0];
  }

  static std::vector<std::array<double,4>> CalculateHDynamicsPolymorphic(const Norm& res, const Norm& mut, double nu, const std::array<double,4>& h_init = {0.5,0.5,0.5,0.5}, size_t output_interval = 0, double dt = 0.01, size_t n_iter = 10000) {
    auto [rr_GG, rr_GB, rr_BG, rr_BB] = RP(res.P, res.R);  // donor: resident, observer: resident
    auto [rm_GG, rm_GB, rm_BG, rm_BB] = RP(res.P, mut.R);  // donor: resident, observer: mutant
    auto [mr_GG, mr_GB, mr_BG, mr_BB] = RP(mut.P, res.R);  // donor: mutant, observer: resident
    auto [mm_GG, mm_GB, mm_BG, mm_BB] = RP(mut.P, mut.R);  // donor: mutant, observer: mutant

    auto hdot = [=](const std::vector<double>& hv, std::vector<double>& dhv) {
      double h_rr = hv[0], h_rm = hv[1], h_mr = hv[2], h_mm = hv[3];
      double dh_rr, dh_rm, dh_mr, dh_mm;
      dh_rr  = nu         * ( h_rm * h_rm * rr_GG + h_rm * (1.0 - h_rm) * rr_GB + (1.0 - h_rm) * h_rm * rr_BG + (1.0 - h_rm) * (1.0 - h_rm) * rr_BB )  // (resident,mutant,resident)
             + (1.0 - nu) * ( h_rr * h_rr * rr_GG + h_rr * (1.0 - h_rr) * rr_GB + (1.0 - h_rr) * h_rr * rr_BG + (1.0 - h_rr) * (1.0 - h_rr) * rr_BB )  // (resident,resident,resident)
             - h_rr;
      dh_rm  = nu         * ( h_mm * h_rm * mr_GG + h_mm * (1.0 - h_rm) * mr_GB + (1.0 - h_mm) * h_rm * mr_BG + (1.0 - h_mm) * (1.0 - h_rm) * mr_BB )  // (mutant,mutant,resident)
             + (1.0 - nu) * ( h_mr * h_rr * mr_GG + h_mr * (1.0 - h_rr) * mr_GB + (1.0 - h_mr) * h_rr * mr_BG + (1.0 - h_mr) * (1.0 - h_rr) * mr_BB )  // (mutant,resident,resident)
             - h_rm;
      dh_mr  = nu         * ( h_rm * h_mm * rm_GG + h_rm * (1.0 - h_mm) * rm_GB + (1.0 - h_rm) * h_mm * rm_BG + (1.0 - h_rm) * (1.0 - h_mm) * rm_BB )  // (resident,mutant,mutant)
             + (1.0 - nu) * ( h_rr * h_mr * rm_GG + h_rr * (1.0 - h_mr) * rm_GB + (1.0 - h_rr) * h_mr * rm_BG + (1.0 - h_rr) * (1.0 - h_mr) * rm_BB )  // (resident,resident,mutant)
             - h_mr;
      dh_mm  = nu         * ( h_mm * h_mm * mm_GG + h_mm * (1.0 - h_mm) * mm_GB + (1.0 - h_mm) * h_mm * mm_BG + (1.0 - h_mm) * (1.0 - h_mm) * mm_BB )  // (mutant,mutant,mutant)
             + (1.0 - nu) * ( h_mr * h_mr * mm_GG + h_mr * (1.0 - h_mr) * mm_GB + (1.0 - h_mr) * h_mr * mm_BG + (1.0 - h_mr) * (1.0 - h_mr) * mm_BB )  // (mutant,resident,mutant)
             - h_mm;
      dhv[0] = dh_rr;
      dhv[1] = dh_rm;
      dhv[2] = dh_mr;
      dhv[3] = dh_mm;
    };

    std::vector<double> hvec(h_init.begin(), h_init.end());
    if (output_interval == 0) {
      output_interval = n_iter;
      hvec = SolveByRungeKutta(hdot, hvec, dt, output_interval);
      return {std::array<double,4>{hvec[0], hvec[1], hvec[2], hvec[3]}};
    }
    else {
      std::vector<std::array<double,4>> ans;
      for (size_t t = 0; t < n_iter; t+=output_interval) {
        hvec = SolveByRungeKutta(hdot, hvec, dt, output_interval);
        std::cerr << t << ' ' << hvec[0] << ' ' << hvec[1] << ' ' << hvec[2] << ' ' << hvec[3] << std::endl;
        ans.emplace_back(std::array<double,4>{hvec[0], hvec[1], hvec[2], hvec[3]});
      }
      return ans;
    }
  }

  static std::array<double,4> CooperationLevelsPolymorphic(const Norm& res, const Norm& mut, double nu, std::optional<std::array<double,4>> h = std::nullopt) {
    // res: resident norm
    // mut: mutant norm
    // nu: fraction of mutants (0, 1)
    // return: [pc_rr, pc_rm, pc_mr, pc_mm]
    if (!h.has_value()) {
      h = HscorePolymorphic(res, mut, nu);
    }
    constexpr Reputation B = Reputation::B, G = Reputation::G;
    double pc_rr = h.value()[0] * res.CProb(G) + (1.0 - h.value()[0]) * res.CProb(B);
    double pc_rm = h.value()[1] * res.CProb(G) + (1.0 - h.value()[1]) * res.CProb(B);
    double pc_mr = h.value()[2] * mut.CProb(G) + (1.0 - h.value()[2]) * mut.CProb(B);
    double pc_mm = h.value()[3] * mut.CProb(G) + (1.0 - h.value()[3]) * mut.CProb(B);
    std::array<double,4> pcs = {pc_rr, pc_rm, pc_mr, pc_mm};
    return pcs;
  }

  static std::array<double,2> Payoffs(const Norm& res, const Norm& mut, double nu, double benefit, std::optional<std::array<double,4>> pcs = std::nullopt) {
    // res: resident norm
    // mut: mutant norm
    // nu: fraction of mutants (0, 1)
    // benefit: benefit of cooperation
    // pcs: cooperation levels [pc_rr, pc_rm, pc_mr, pc_mm] (optional)
    // return: [payoff_res, payoff_mut]
    if (!pcs.has_value()) {
      pcs = CooperationLevelsPolymorphic(res, mut, nu);
    }
    auto [pc_rr, pc_rm, pc_mr, pc_mm] = pcs.value();
    double payoff_res = benefit * ( nu * pc_mr + (1.0 - nu) * pc_rr ) - ( nu * pc_rm + (1.0 - nu) * pc_rr );
    double payoff_mut = benefit * ( nu * pc_mm + (1.0 - nu) * pc_rm ) - ( nu * pc_mm + (1.0 - nu) * pc_mr );
    return {payoff_res, payoff_mut};
  }

  using vd_t = std::vector<double>;
  static vd_t SolveByRungeKutta(const std::function<void(const vd_t&,vd_t&)>& func, const vd_t& init, double dt, size_t n_iter) {
    const size_t N = init.size();
    vd_t ht = init;
    vd_t k1(N, 0.0), arg2(N, 0.0), k2(N, 0.0), arg3(N, 0.0), k3(N, 0.0), arg4(N, 0.0), k4(N, 0.0);
    for (size_t t = 0; t < n_iter; t++) {
      func(ht, k1);
      for(int i = 0; i < N; i++) {
        k1[i] *= dt;
        arg2[i] = ht[i] + 0.5 * k1[i];
      }
      func(arg2, k2);
      for(int i = 0; i < N; i++) {
        k2[i] *= dt;
        arg3[i] = ht[i] + 0.5 * k2[i];
      }
      func(arg3, k3);
      for(int i = 0; i < N; i++) {
        k3[i] *= dt;
        arg4[i] = ht[i] + k3[i];
      }
      func(arg4, k4);
      for(int i = 0; i < N; i++) {
        k4[i] *= dt;
      }
      for (int i = 0; i < N; i++) {
        ht[i] += (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) / 6.0;
      }
      // std::cerr << t << ' ' << ht[0] << ' ' << ht[1] << ' ' << ht[2] << ' ' << ht[3] << std::endl;
    }
    return ht;
  }
};

#endif // SOLITARY_OBS_GAME_HPP
