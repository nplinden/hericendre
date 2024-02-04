#ifndef NFY_HPP_INCLUDED
#define NFY_HPP_INCLUDED
#include <pugixml.hpp>
#include <vector>
#include <map>
#include <string>

class NeutronFissionYield{
    public:
        NeutronFissionYield(const pugi::xml_node& NFYNode) ;
        std::vector<double> energies_ ;
        std::map<std::string, std::map<std::string, double>> fy_ ;
} ;

#endif