#include "results.h"
#include <fmt/os.h>
#include <highfive/H5Easy.hpp>

Results::Results() = default;

Results::Results(const std::vector<std::vector<double>> &cc,
                 const std::vector<std::string> &nuclides, const std::vector<double> &times)
{
    cc_ = cc;
    nuclides_ = nuclides;
    times_ = times;
};

Results::Results(const std::vector<Eigen::VectorXd> &cc,
                 const std::vector<std::string> &nuclides, const std::vector<double> &times)
{
    for (auto &vec : cc)
    {
        cc_.emplace_back(vec.data(), vec.data() + vec.size());
    }
    nuclides_ = nuclides;
    times_ = times;
};

void Results::to_csv(const std::string &path, const bool ignore_zeros) const
{
    auto out = fmt::output_file(path);

    // HEADER
    for (const auto &t : this->times_)
        out.print(",{}", t);
    out.print("\n");

    // CONTENT
    for (size_t inuc = 0; inuc < this->nuclides_.size(); inuc++)
    {
        std::string line = fmt::format("{:s}", this->nuclides_[inuc]);
        bool zero_flag = true;
        for (const auto &concentrations : this->cc_)
        {
            line += fmt::format(",{:e}", concentrations[inuc]);
            if (concentrations[inuc] != 0.)
                zero_flag = false;
        }
        line += fmt::format("\n");
        if (!zero_flag || !ignore_zeros)
        {
            out.print("{}", line);
        }
    }
};

void Results::to_hdf5(H5Easy::File &file) const
{
    H5Easy::dump(file, "CONCENTRATIONS", this->cc_);
    H5Easy::dump(file, "TIMES", this->times_);
    H5Easy::dump(file, "NUCLIDES", this->nuclides_);
}
