#ifndef CHAIN_HPP_INCLUDED
#define CHAIN_HPP_INCLUDED
#include <string>
#include <vector>
#include <nuclide.hpp>
#include <decay.hpp>
#include <nreaction.hpp>
#include <Eigen/Sparse>

using DecayPtr = std::shared_ptr<Decay>;
using NReactionPtr = std::shared_ptr<NReaction>;
using NuclidePtr = std::shared_ptr<Nuclide>;
class Chain{
    public:
        Chain(const char* path) ;
        std::vector<NuclidePtr> nuclides_ ;
        std::vector<DecayPtr> decays_ ;
        std::vector<NReactionPtr> nreactions_ ;
        NuclidePtr find(int nucid) const;
        NuclidePtr find(std::string name) const;
        int nuclide_index(int nucid) const;
        int nuclide_index(std::string nuclide_name) const;
        Eigen::SparseMatrix<double> decayMatrix() const ;
} ;

#endif