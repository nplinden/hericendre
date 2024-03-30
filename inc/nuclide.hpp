#ifndef NUCLIDE_HPP_INCLUDED
#define NUCLIDE_HPP_INCLUDED
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <pugixml.hpp>
#include <fmt/core.h>
#include "source.hpp"
// #include "nfy.hpp"

class Decay;
class NReaction;
class Fission;
using SourcePtr = std::shared_ptr<Source>;
// using NFYPtr = std::shared_ptr<NeutronFissionYield>;
class Nuclide{
    public:
        // Constructors
        Nuclide(std::string name, double dconst) ;
        Nuclide(const pugi::xml_node& nuclideNode) ;

        // Members functions
        void addNode(pugi::xml_node& rootnode) ;
        bool operator<(const Nuclide& other) const ;
        std::string str() ;
        static std::tuple<int, int, int> getZam(std::string name) ;

        // Member variables
        // NFYPtr neutronFissionYields_ ;
        std::vector<std::shared_ptr<Decay>> decays_ ;
        std::vector<std::shared_ptr<Decay>> decaysUp_ ;
        std::vector<std::shared_ptr<NReaction>> reactions_ ;
        std::vector<std::shared_ptr<NReaction>> reactionsUp_ ;
        std::vector<SourcePtr> sources_ ;
        double dconst_ ;
        double denergy_ ;
        std::string name_ ;
        std::string nfyParent_ = "" ;
        int idInChain ;
        int nreac_ ;
        int ndecay_ ;
        static const std::map<std::string, int> ELEMENTS ;
        int z_ ;
        int a_ ;
        int m_ ;
        int zam_ ;
} ;

template <>
class fmt::formatter<Nuclide> {
public:
  constexpr auto parse (format_parse_context& ctx) { return ctx.begin(); }
  template <typename Context>
  constexpr auto format (Nuclide const& nuc, Context& ctx) const {
      return format_to(ctx.out(), "{}", nuc.name_);
  }
};
#endif
