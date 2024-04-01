#include "pugixml.hpp"
#include <Eigen/Sparse>
#include <chain.h>
#include <decay.h>
#include <fmt/core.h>
#include <fmt/os.h>
#include <fmt/ranges.h>
#include <iostream>
#include <map>
#include <nuclide.h>
#include <numeric>
#include <solver.h>
#include <vector>

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

enum class RunMode { PROCESS, DEPLETE, BURN };

Eigen::VectorXd deplete(Chain &chain) {
  std::map<std::string, double> ccMap;

  Eigen::VectorXd N(chain.nuclides_.size());
  for (int i = 0; i < chain.nuclides_.size(); i++)
    N(i) = 1.;
  Solver solver;
  return solver.run(chain, N, 86400);
}

int main(int argc, char *argv[]) {
  Chain chain("/home/nlinden/workspace/hericendre/data/chain_endfb71_sfr.xml");
  // Chain
  // chain("/home/nlinden/workspace/hericendre/data/chain-jeff33-long.xml");
  auto nuclides = chain.reachable("Pu239");

  std::map<std::string, double> ccMap;
  ccMap["Pu239"] = 1.e24;

  Solver solver;
  // N = solver.run(chain, N, 86400);
  std::vector<double> times = {
      1e10, 1e11, 1e12, 1e13, 1e14,
  };
  auto N = solver.run(chain, ccMap, times);

  // for (int i = 0; i < N.size(); i++) {
  //   if (N[i] != 0)
  //     fmt::print("{:8s},{:.4e},{: .4e}\n", chain.nuclides_[i]->name_,
  //                chain.nuclides_[i]->dconst_, N[i]);
  // }
  // double sum = std::accumulate(N.begin(), N.end(), 0.);
  // fmt::print("sum = {:.4e}\n", sum);

  return 0;
}
