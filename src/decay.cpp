#include "decay.hpp"
#include "nuclide.hpp"


Decay::Decay(const pugi::xml_node& decayNode, NuclidePtr parent){
    type_ = decayNode.attribute("type").value() ;
    targetName_ = decayNode.attribute("target").value() ;

    std::string br = decayNode.attribute("branching_ratio").value() ;
    if (!br.empty())
        branchingRatio_ = stod(br) ;
    else
        branchingRatio_ = 1. ;

    parent_ = parent ;
} ;