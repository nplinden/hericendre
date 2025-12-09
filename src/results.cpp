#include "results.h"
#include <fmt/os.h>
#include <highfive/H5Easy.hpp>

Results::Results() = default;

Results::Results(const std::vector<std::vector<double>> &cc,
                 const std::vector<std::string> &nuclides,
                 const std::vector<double> &times)
    : cc_(cc), nuclides_(nuclides), times_(times){}

Results::Results(const std::vector<Eigen::VectorXd> &cc,
                 const std::vector<std::string> &nuclides,
                 const std::vector<double> &times): nuclides_(nuclides), times_(times) {
  cc_.reserve(cc.size());
  for (const auto &vec : cc) {
    cc_.emplace_back(vec.data(), vec.data() + vec.size());
  }
}

void Results::to_csv(const std::string &path, bool ignore_zeros) const {
  auto out = fmt::output_file(path);

  // HEADER
  out.print("nuclide");
  for (const auto &t : times_) {
    out.print(",{}", t);
  }
  out.print("\n");

  // CONTENT
  for (size_t inuc = 0; inuc < nuclides_.size(); inuc++) {
    bool has_nonzero = false;
    
    // Check if this nuclide has any non-zero values
    if (ignore_zeros) {
      for (const auto &concentrations : cc_) {
        if (concentrations[inuc] != 0.0) {
          has_nonzero = true;
          break;
        }
      }
      if (!has_nonzero) {
        continue;
      }
    }
    
    // Write the row
    out.print("{}", nuclides_[inuc]);
    for (const auto &concentrations : cc_) {
      out.print(",{:e}", concentrations[inuc]);
    }
    out.print("\n");
  }
}


void Results::to_hdf5(H5Easy::File &file) const {
  H5Easy::dump(file, "CONCENTRATIONS", this->cc_);
  H5Easy::dump(file, "TIMES", this->times_);
  H5Easy::dump(file, "NUCLIDES", this->nuclides_);
}
