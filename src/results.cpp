#include "results.h"
#include <fmt/os.h>

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

void Results::to_csv(std::string path, bool ignore_zeros){
  
  auto out = fmt::output_file(path);

  // HEADER
  for (const auto& t: this->times_)
    out.print(",{}", t);
  out.print("\n");

  // CONTENT
  bool zero_flag;
  std::string line;
  for (size_t inuc = 0; inuc < this->nuclides_.size(); inuc++){
    line = fmt::format("{:s}", this->nuclides_[inuc]);
    zero_flag = true;
    for (const auto& concentrations: this->cc_){
      line += fmt::format(",{:e}", concentrations[inuc]);
      if (concentrations[inuc] != 0.)
        zero_flag = false;
    }
    line += fmt::format("\n");
    if (!zero_flag || !ignore_zeros){
      out.print("{}", line);
    }
  }
};