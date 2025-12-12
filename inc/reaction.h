#ifndef REACTION_HPP_INCLUDED
#define REACTION_HPP_INCLUDED
#include <memory>
#include <pugixml.hpp>

class Nuclide;
using NuclidePtr = std::shared_ptr<Nuclide>;

class Reaction
{
public:
    Reaction(const pugi::xml_node &reactionNode, const NuclidePtr &parent);
    NuclidePtr parent_;
    NuclidePtr target_;
    std::string targetName_;
    std::string type_;

    double xs_;
    double branchingRatio_;

private:
    double Q_;
};
#endif // REACTION_HPP_INCLUDED