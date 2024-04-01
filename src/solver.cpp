#include <Eigen/SparseLU>
#include <fmt/core.h>
#include <solver.h>

Solver::Solver() {}

Eigen::VectorXd Solver::run(const Chain &chain, Eigen::VectorXd ccVector,
                            double dt, double cutoff) {
  int n_nuclides = chain.nuclides_.size();
  SpComplex M = chain.decayMatrix().cast<cdouble>();

  std::vector<TrComplex> t_triplets;
  t_triplets.reserve(n_nuclides);
  for (int j = 0; j < n_nuclides; j++)
    t_triplets.push_back(TrComplex(j, j, std::complex(1., 0.)));
  SpComplex Identity(n_nuclides, n_nuclides);
  Identity.setFromTriplets(t_triplets.begin(), t_triplets.end());

  Eigen::VectorX<cdouble> N = ccVector.cast<cdouble>();
  Eigen::SparseLU<SpComplex, Eigen::COLAMDOrdering<int>> solver;
  for (int i = 0; i < theta48.size(); i++) {
    cdouble theta = theta48[i];
    cdouble alpha = alpha48[i];

    SpComplex A = M * dt - theta * Identity;
    solver.analyzePattern(A);
    solver.factorize(A);
    Eigen::VectorX<cdouble> x = (alpha * solver.solve(N));
    N += 2 * x.real();
  }

  Eigen::VectorX<double> realN = alpha48_0 * N.real();

  for (int i = 0; i < realN.size(); i++) {
    if (realN[i] < cutoff)
      realN[i] = 0.;
  }

  return realN;
}

Eigen::VectorXd Solver::run(const Chain &chain,
                            std::map<std::string, double> ccMap, double dt,
                            double cutoff) {
  Eigen::VectorXd N(chain.nuclides_.size());
  for (int i; i < chain.nuclides_.size(); i++)
    N(i) = 0.;

  for (auto const &[key, val] : ccMap) {
    int inuc = chain.nuclide_index(key);
    N(inuc) = val;
  }
  return run(chain, N, dt, cutoff);
}

std::vector<Eigen::VectorXd> Solver::run(const Chain &chain,
                                         Eigen::VectorXd ccVector,
                                         std::vector<double> times,
                                         double cutoff) {
  std::vector<double> dts;
  dts.push_back(times[0]);
  for (int it = 1; it < dts.size(); it++)
    dts.push_back(times[it] - times[it - 1]);

  std::vector<Eigen::VectorXd> results;
  Eigen::VectorXd N(ccVector);

  results.push_back(N);
  for (const auto &dt : dts) {
    N = run(chain, N, dt, cutoff);
    results.push_back(N);
  }
  return results;
}

std::vector<Eigen::VectorXd> Solver::run(const Chain &chain,
                                         std::map<std::string, double> ccMap,
                                         std::vector<double> times,
                                         double cutoff) {
  Eigen::VectorXd N(chain.nuclides_.size());
  for (int i; i < chain.nuclides_.size(); i++)
    N(i) = 0.;

  for (auto const &[key, val] : ccMap) {
    int inuc = chain.nuclide_index(key);
    N(inuc) = val;
  }

  return run(chain, N, times, cutoff);
}
