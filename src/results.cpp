#include "results.h"

Results::Results(){};
Results::Results(std::vector<std::vector<double>> cc,
                 std::vector<std::string> nuclides, std::vector<double> times) {
  cc_ = cc;
  nuclides_ = nuclides;
  times_ = times;
};
Results::Results(std::vector<Eigen::VectorXd> cc,
                 std::vector<std::string> nuclides, std::vector<double> times) {
  for (auto &vec : cc) {
    cc_.push_back(std::vector<double>(vec.data(), vec.data() + vec.size()));
  }
  nuclides_ = nuclides;
  times_ = times;
};
