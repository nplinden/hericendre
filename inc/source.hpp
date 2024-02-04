#ifndef SOURCE_HPP_INCLUDED
#define SOURCE_HPP_INCLUDED
#include <string>
#include <vector>
#include <pugixml.hpp>

class Source{
    public:
        Source(const pugi::xml_node& sourceNode) ;
        std::string type_ ;
        std::string particle_ ;
        std::string interpolation_ ;
        std::vector<double> energy_ ;
        std::vector<double> intensity_ ;
        std::vector<Source> pairs_ ;
        std::vector<double> pair_probabilities_ ;
} ;

#endif