#ifndef DECAYSOLVER_HPP_INCLUDED
#define DECAYSOLVER_HPP_INCLUDED
#include "chain.h"
#include "results.h"
#include <map>
#include <string>
#include <highfive/H5Easy.hpp>

class DecaySolver {
public:
    DecaySolver();

    void compute_coeffs(Chain &chain, std::map<std::string, double> ccMap);

    std::vector<std::vector<double> > run(Chain &chain,
                                          const std::map<std::string, double> &ccMap,
                                          std::vector<double> times);

    void to_hdf5(H5Easy::File &file) const;

    Chain chain_;
    std::map<std::string, double> ccMap_;

    std::map<size_t, double> Ns;
    std::map<size_t, std::map<size_t, double> > F;
    Results results_;
};

#endif
