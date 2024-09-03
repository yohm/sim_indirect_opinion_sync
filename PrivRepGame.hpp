#ifndef PRIVATE_REP_GAME_HPP
#define PRIVATE_REP_GAME_HPP

#include <iostream>
#include <vector>
#include <random>
#include <cassert>
#include "Norm.hpp"

class PrivateRepGame {
public:
  using population_t = std::vector<std::pair<Norm, size_t>>;  // vector of NormID & its size

  PrivateRepGame(population_t pop, uint64_t _seed) : population(std::move(pop)), rnd(_seed), uni(0.0, 1.0) {
    for (size_t i = 0; i < population.size(); i++) {
      auto kv = population[i];
      for (size_t ii = 0; ii < kv.second; ii++) {
        norms.emplace_back(kv.first);
        strategy_idx.emplace_back(i);
      }
    }
    N = norms.size();
    M.assign(N*N, Reputation::G);

    // randomize initial image
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < N; j++) {
        M[i*N+j] = rnd() % 2 == 0 ? Reputation::G : Reputation::B;
      }
    }

    // initialize coop_count with NxN matrix of 0
    coop_count.assign(N*N, 0);
    game_count.assign(N*N, 0);
    good_count.assign(N*N, 0);
    good_count_times = 0;
    gi_sum.assign(N, 0.0);
    gi_squared_sum.assign(N, 0.0);
  }

  static std::pair<std::vector<size_t>,std::vector<size_t>> RandomNonIdenticalPermutations(size_t N, std::mt19937_64& rnd) {
    // construct two random permutations of the numbers ranging from 0 to N-1,
    // while ensuring that the i-th number in the first permutation is never identical to the i-th number in the second permutation
    std::vector<size_t> perm1(N), perm2(N);
    for (size_t i = 0; i < N; i++) {
      perm1[i] = i;
      perm2[i] = i;
    }
    std::shuffle(perm1.begin(), perm1.end(), rnd);
    std::shuffle(perm2.begin(), perm2.end(), rnd);

    // perm2[i] should be different from perm1[i] for any i
    // if not, swap perm2[i] with perm2[j] where j is the smallest index such that perm2[j] != perm1[i]
    for (size_t i = 0; i < N; i++) {
      if (perm1[i] == perm2[i]) {
        size_t j = (i+1)%N;
        while (perm1[j] == perm2[i] || perm1[i] == perm2[j]) {
          j++;
          j %= N;
        }
        std::swap(perm2[i], perm2[j]);
      }
    }

    return std::make_pair(perm1, perm2);
  }

  // t_max : number of steps
  // q : observation probability
  void UpdateGossiping(size_t t_max, double q, double mu_impl, double mu_assess, double tau, bool count_good) {
    for (size_t t = 0; t < t_max; t++) {

      size_t NG = 1 * N;
      std::vector<Action> actions(NG);
      std::vector< std::array<size_t,3> > pairs(NG);
      for (size_t tt = 0; tt < NG; tt++) {
        std::uniform_int_distribution<size_t> uni_n(0, N-1);
        size_t donor = uni_n(rnd);
        std::uniform_int_distribution<size_t> uni_n1(1, N-1);
        size_t recip = (donor + uni_n1(rnd)) % N;
        size_t obs = (donor + uni_n1(rnd)) % N;
        assert(donor != recip);

        Action A = DecideAction(donor, recip, mu_impl);
        actions[tt] = A;
        pairs[tt] = {donor, recip, obs};

        if (A == Action::C) {
          coop_count[donor * N + recip]++;
        }
        game_count[donor * N + recip]++;
      }

      std::vector<Reputation> assessments(NG);
      for (size_t tt = 0; tt < NG; tt++) {
        auto [donor,recip,obs] = pairs[tt];
        // randomly choose an observer excluding the donor
        assessments[tt] = AssessmentByObserverAboutDonor(obs, donor, recip, actions[tt], mu_assess);
      }
      for (size_t tt = 0; tt < NG; tt++) {
        auto [donor,recip,obs] = pairs[tt];
        M[obs * N + donor] = assessments[tt];
      }

      // if (t % N == N-1) {
      //   auto num_gossip_events = static_cast<size_t>( std::round(0.5 * tau * N * NG * N));
      //   ExchangeOpinionsRandom(num_gossip_events);
      // }
      auto num_gossip_events = static_cast<size_t>( std::round(0.5 * tau * N * NG));
      ExchangeOpinionsRandom(num_gossip_events);

      // count good
      if (count_good) {
        for (size_t i = 0; i < N; i++) {
          int good_i = 0;
          for (size_t j = 0; j < N; j++) {
            if (M[j * N + i] == Reputation::G) {
              good_count[j * N + i]++;
              if (j != i) {  // exclude the diagonal elements when calculating goodness
                good_i++;
              }
            }
          }
          double gi = static_cast<double>(good_i) / static_cast<double>(N - 1);
          gi_sum[i] += gi;
          gi_squared_sum[i] += gi * gi;
        }
        good_count_times++;
      }
    }
  }


  // t_max : number of steps
  // q : observation probability
  void Update(size_t t_max, double q, double mu_impl, double mu_assess, double tau, bool count_good) {
    for (size_t t = 0; t < t_max; t++) {
      auto [perm1, perm2] = RandomNonIdenticalPermutations(N, rnd);

      for (size_t tt=0; tt < N; tt++) {
        // randomly choose donor & recipient
        size_t donor = perm1[tt];
        size_t recip = perm2[tt];
        assert(donor != recip);

        Action A = DecideAction(donor, recip, mu_impl);

        if (A == Action::C) {
          coop_count[donor * N + recip]++;
        }
        game_count[donor * N + recip]++;

        // action observation phase
        if (q > 0.0) {
          for (size_t obs = 0; obs < N; obs++) {
            if (obs == donor || q == 1.0 || R01() < q) {  // observe with probability q
              M[obs * N + donor] = AssessmentByObserverAboutDonor(obs, donor, recip, A, mu_assess);
            }
          }
        }
        else {  // solitary observation
          // randomly choose an observer excluding the donor
          std::uniform_int_distribution<size_t> uni_n(1, N-1);
          size_t obs = (donor + uni_n(rnd)) % N;
          M[obs * N + donor] = AssessmentByObserverAboutDonor(obs, donor, recip, A, mu_assess);

          // self-image of the donor is also updated
          M[donor * N + donor] = AssessmentByObserverAboutDonor(donor, donor, recip, A, mu_assess);
        }

        // exchange opinions via gossip
        // auto num_gossip_events = static_cast<size_t>( std::round(0.5 * tau * N) );
        // ExchangeOpinionsRandom(num_gossip_events);
        // ExchangeOpinions(donor, num_gossip_events);
      }

      //if (t % N == N-1) {
        auto num_gossip_events = static_cast<size_t>( std::round(0.5 * tau * N * N) );
        ExchangeOpinionsRandom(num_gossip_events);
      //}

      // count good
      if (count_good) {
        for (size_t i = 0; i < N; i++) {
          int good_i = 0;
          for (size_t j = 0; j < N; j++) {
            if (M[j * N + i] == Reputation::G) {
              good_count[j * N + i]++;
              if (j != i) {  // exclude the diagonal elements when calculating goodness
                good_i++;
              }
            }
          }
          double gi = static_cast<double>(good_i) / static_cast<double>(N-1);
          gi_sum[i] += gi;
          gi_squared_sum[i] += gi * gi;
        }
        good_count_times++;
      }
    }
  }

private:
  Action DecideAction(size_t donor, size_t recip, double mu_impl) {
    double c_prob = norms[donor].CProb(M[donor * N + recip]);
    Action A;
    if (c_prob == 1.0 && mu_impl == 0.0) {
      A = Action::C;
    } else if (c_prob == 0.0) {
      A = Action::D;
    } else {
      A = (R01() < c_prob * (1.0-mu_impl)) ? Action::C : Action::D;
    }
    return A;
  }

  Reputation AssessmentByObserverAboutDonor(size_t obs, size_t donor, size_t recip, Action A, double mu_assess) {
    double g_prob = norms[obs].R.GProb(M[obs * N + recip], A);
    Reputation rep;
    if (g_prob == 1.0) {
      rep = Reputation::G;
    } else if (g_prob == 0.0) {
      rep = Reputation::B;
    } else {
      rep = (R01() < g_prob) ? Reputation::G : Reputation::B;
    }
    // assessment error
    if (mu_assess > 0.0 && R01() < mu_assess) {
      rep = FlipReputation(rep);
    }
    return rep;
  }

  void ExchangeOpinionsRandom(size_t T) {
    // exchange opinions about i T times
    for (size_t t = 0; t < T; t++) {
      // draw three unique random numbers i,j,k
      std::uniform_int_distribution<size_t> uni_n(0, N - 1);
      size_t i = uni_n(rnd);
      size_t j = uni_n(rnd);
      size_t k = uni_n(rnd);
      while (i == j || i == k || j == k) {
        i = uni_n(rnd);
        j = uni_n(rnd);
        k = uni_n(rnd);
      }
      // copy m[i][j] to copy[k][j]
      M[k * N + j] = M[i * N + j];
    }
  }

  void ExchangeOpinions(size_t i, size_t T) {
    // find the fraction goodness of i
    size_t num_good = 0;
    for (size_t j = 0; j < N; j++) {
      if (j == i) continue;
      if (M[j * N + i] == Reputation::G) {
        num_good++;
      }
    }
    // exchange opinions about i
    std::uniform_int_distribution<size_t> uni_n(1, N-1);
    for (size_t t = 0; t < T; t++) {
      size_t j = (i + uni_n(rnd)) % N;
      Reputation org = M[j * N + i];
      // we want to find the fraction of good reputations excluding i and j
      if (org == Reputation::G) {
        num_good--;
      }
      double goodness = static_cast<double>(num_good) / static_cast<double>(N-2);
      M[j * N + i] = (R01() < goodness) ? Reputation::G : Reputation::B;
      if (M[j * N + i] == Reputation::G) {
        num_good++;
      }
    }
  }

public:
  using count_t = std::vector<size_t>;

  // number of C between i-donor and j-recipient
  const count_t& CoopCount() const { return coop_count; }

  // number of games between i-donor and j-recipient
  const count_t& GameCount() const { return game_count; }

  // system-wide cooperation level
  double SystemWideCooperationLevel() const {
    size_t coop = 0;
    size_t total = 0;
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < N; j++) {
        coop += coop_count[i*N+j];
        total += game_count[i*N+j];
      }
    }
    return static_cast<double>(coop) / static_cast<double>(total);
  }

  // donation levels of individuals
  // <probability of receiving benefit, probability of paying cost>
  // payoff
  std::vector<std::pair<double,double>> IndividualCooperationLevels() const {
    std::vector<std::pair<double,double>> coop_rates(N, {0.0, 0.0});
    for (size_t i = 0; i < N; i++) {
      size_t coop_total_out = 0, game_total_out = 0;
      for (size_t j = 0; j < N; j++) {
        if (i == j) {
          continue;
        }
        coop_total_out += coop_count[i*N+j];
        game_total_out += game_count[i*N+j];
      }
      size_t coop_total_in = 0.0, game_total_in = 0.0;
      for (size_t j = 0; j < N; j++) {
        if (i == j) {
          continue;
        }
        coop_total_in += coop_count[j*N+i];
        game_total_in += game_count[j*N+i];
      }
      double in = static_cast<double>(coop_total_in) / static_cast<double>(game_total_in);
      double out = static_cast<double>(coop_total_out) / static_cast<double>(game_total_out);
      coop_rates[i] = std::make_pair(in, out);
    }
    return coop_rates;
  }

  // norm-wise cooperation levels
  // return: vector c_levels
  //   c_levels[i][j] : cooperation level of i-th norm toward j-th norm
  std::vector<std::vector<double>> NormCooperationLevels() const {
    size_t n_norms = population.size();
    std::vector<std::vector<size_t>> coop_by_norm(n_norms, std::vector<size_t>(n_norms, 0) );
    std::vector<std::vector<size_t>> total_by_norm(n_norms, std::vector<size_t>(n_norms, 0) );
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < N; j++) {
        size_t i_norm = strategy_idx[i];
        size_t j_norm = strategy_idx[j];
        coop_by_norm[i_norm][j_norm] += coop_count[i*N+j];
        total_by_norm[i_norm][j_norm] += game_count[i*N+j];
      }
    }
    std::vector<std::vector<double>> c_levels(n_norms, std::vector<double>(n_norms, 0.0));
    for (size_t i_norm = 0; i_norm < n_norms; i_norm++) {
      for (size_t j_norm = 0; j_norm < n_norms; j_norm++) {
        c_levels[i_norm][j_norm] = static_cast<double>(coop_by_norm[i_norm][j_norm]) / static_cast<double>(total_by_norm[i_norm][j_norm]);
      }
    }
    return c_levels;
  }

  // norm-wise good count
  // return: vector average reputation
  //   c_levels[i][j] : average reputation of j-th norm from the viewpoint of i-th norm
  std::vector<std::vector<double>> NormAverageReputation() const {
    size_t n_norms = population.size();
    std::vector<std::vector<double>> avg_rep(n_norms, std::vector<double>(n_norms, 0.0));

    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < N; j++) {
        size_t i_norm = strategy_idx[i];
        size_t j_norm = strategy_idx[j];
        avg_rep[i_norm][j_norm] += static_cast<double>(good_count[i*N+j]);
      }
    }

    for (size_t i_norm = 0; i_norm < n_norms; i_norm++) {
      for (size_t j_norm = 0; j_norm < n_norms; j_norm++) {
        size_t c = population[i_norm].second * population[j_norm].second * good_count_times;
        avg_rep[i_norm][j_norm] /= static_cast<double>(c);
      }
    }

    return avg_rep;
  }

  // norm-wise goodness count
  // return: vector of pairs {<gi>, <gi^2>}
  std::vector<std::pair<double,double>> NormAverageGoodness() const {
    size_t n_norms = population.size();
    std::vector<std::pair<double,double>> avg_goodness(n_norms, {0.0, 0.0});
    for (size_t i = 0; i < N; i++) {
      size_t i_norm = strategy_idx[i];
      avg_goodness[i_norm].first += gi_sum[i];
      avg_goodness[i_norm].second += gi_squared_sum[i];
    }
    for (size_t i_norm = 0; i_norm < n_norms; i_norm++) {
      double denom = population[i_norm].second * good_count_times;
      avg_goodness[i_norm].first /= denom;
      avg_goodness[i_norm].second /= denom;
    }
    return avg_goodness;
  }

  void ResetCounts() {
    // reset coop_count and fill with 0
    std::fill(coop_count.begin(), coop_count.end(), 0);
    std::fill(game_count.begin(), game_count.end(), 0);
    std::fill(good_count.begin(), good_count.end(), 0);
    good_count_times = 0;
  }

  // print image matrix
  void PrintImage(std::ostream& out) const {
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < N; j++) {
        out << ((M[i*N+j]==Reputation::G)?'.':'x');
      }
      out << "\n";
    }
  }

  population_t Population() const {
    return population;
  }

private:
  const population_t population;
  size_t N;
  std::mt19937_64 rnd;
  std::uniform_real_distribution<double> uni;
  std::vector<Norm> norms;
  std::vector<size_t> strategy_idx;  // [0,0,...,0, 1,1,....,1, 2,2,......,2]
  std::vector<Reputation> M;
  count_t coop_count;  // number of C between i-donor and j-recipient
  count_t game_count;  // number of games between i-donor and j-recipient
  count_t good_count;  // number of good reputations
  size_t good_count_times; // number of time steps for which good_count is counted
  std::vector<double> gi_sum;  // sum of gi
  std::vector<double> gi_squared_sum;  // sum of gi^2
  // <gi> = gi_sum[i] / good_count_times, <gi^2> = gi_squared_sum[i] / good_count_times
  double R01() { return uni(rnd); }
};

#endif // PRIVATE_REP_GAME_HPP