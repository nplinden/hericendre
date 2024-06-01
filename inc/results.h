#ifndef RESULTS_HPP_INCLUDED
#define RESULTS_HPP_INCLUDED
#include <Eigen/Sparse>
#include <vector>

class Results {
public:
  std::vector<std::string> nuclides_;
  std::vector<double> times_;
  std::vector<std::vector<double>> cc_;

  Results();
  Results(const std::vector<std::vector<double>> &cc,
          const std::vector<std::string> &nuclides, const std::vector<double> &times);
  Results(const std::vector<Eigen::VectorXd> &cc, const std::vector<std::string> &nuclides,
          const std::vector<double> &times);

  void to_csv(const std::string &path, bool ignore_zeros = true);
};

#endif
