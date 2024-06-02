#ifndef DECAYSOLVER_HPP_INCLUDED
#define DECAYSOLVER_HPP_INCLUDED
#include "chain.h"
#include "results.h"
#include <map>
#include <string>

class DecaySolver {
public:
    DecaySolver();

    void compute_coeffs(Chain &chain, std::map<std::string, double> ccMap);

    std::vector<std::vector<double> > run(Chain &chain,
                                          const std::map<std::string, double> &ccMap,
                                          std::vector<double> times);

    std::map<size_t, double> Ns;
    std::map<size_t, std::map<size_t, double> > F;
    Results results_;
};

#endif
