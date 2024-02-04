#ifndef NUCLIDE_HPP_INCLUDED
#define NUCLIDE_HPP_INCLUDED
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <pugixml.hpp>
#include <fmt/core.h>
#include "source.hpp"
#include "nfy.hpp"

class Decay;
class NReaction;
class NeutronFissionYield;
using SourcePtr = std::shared_ptr<Source>;
using NFYPtr = std::shared_ptr<NeutronFissionYield>;
class Nuclide{
    public:
        Nuclide(std::string name, double dconst) ;
        Nuclide(const pugi::xml_node& nuclideNode) ;
        void addNode(pugi::xml_node& rootnode) ;
        std::string str() ;
        static std::tuple<int, int, int> getZam(std::string name) ;
        std::vector<std::shared_ptr<Decay>> decays_ ;
        std::vector<std::shared_ptr<NReaction>> nreactions_ ;
        std::vector<SourcePtr> sources_ ;
        NFYPtr neutronFissionYields_ ;
        bool operator<(const Nuclide& other) const ;
        
        double dconst_ ;
        float denergy_ ;
        std::string name_ ;
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