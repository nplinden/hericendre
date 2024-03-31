#include <Eigen/SparseLU>
#include <fmt/core.h>
#include <solver.hpp>

Solver::Solver() {}

Eigen::VectorXd Solver::run(const Chain &chain,
                            std::map<std::string, double> ccMap, double dt) {
  Eigen::VectorXd N(chain.nuclides_.size());
  for (int i; i < chain.nuclides_.size(); i++)
    N(i) = 0.;

  for (auto const &[key, val] : ccMap) {
    int inuc = chain.nuclide_index(key);
    N(inuc) = val;
  }
  return run(chain, N, dt);
}

Eigen::VectorXd Solver::run(const Chain &chain, Eigen::VectorXd ccVector,
                            double dt) {
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
  return alpha48_0 * N.real();
}
