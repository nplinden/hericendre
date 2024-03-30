#ifndef CHAIN_HPP_INCLUDED
#define CHAIN_HPP_INCLUDED
#include <string>
#include <vector>
#include <nuclide.hpp>
#include <fission.hpp>
#include <decay.hpp>
#include <nreaction.hpp>
#include <Eigen/Sparse>

using DecayPtr = std::shared_ptr<Decay>;
using NReactionPtr = std::shared_ptr<NReaction>;
using NuclidePtr = std::shared_ptr<Nuclide>;
using NFYPtr = std::shared_ptr<Fission>;
class Chain{
    public:
        // Constructors
        Chain(const char* path) ;
        Chain() ;

        // Member Functions
        void write(const char* path) ;
        void restrict(std::vector<std::string> nuclides) ;
        void removeNuclide(std::string nuc) ;
        NuclidePtr find(int nucid) const;
        NuclidePtr find(std::string name) const;
        int nuclide_index(int nucid) const;
        int nuclide_index(std::string name) const;
        Eigen::SparseMatrix<double> decayMatrix() const ;

        // Member variables
        std::vector<NuclidePtr> nuclides_ ;
        std::vector<DecayPtr> decays_ ;
        std::vector<NReactionPtr> reactions_ ;
        std::vector<NFYPtr> fissions_ ;
} ;

#endif