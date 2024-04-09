
#ifndef SOLVER_HPP_INCLUDED
#define SOLVER_HPP_INCLUDED
#include <Eigen/Sparse>
#include <chain.h>
#include <vector>

class Results {
public:
  std::vector<std::string> nuclides_;
  std::vector<double> times_;
  std::vector<std::vector<double>> cc_;

  Results();
  Results(std::vector<std::vector<double>> cc,
          std::vector<std::string> nuclides, std::vector<double> times);
  Results(std::vector<Eigen::VectorXd> cc, std::vector<std::string> nuclides,
          std::vector<double> times);
};

#endif
