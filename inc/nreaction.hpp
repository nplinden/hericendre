#ifndef NREACTION_HPP_INCLUDED
#define NREACTION_HPP_INCLUDED
#include <memory>
#include <pugixml.hpp>
#include <string>

class Nuclide;
using NuclidePtr = std::shared_ptr<Nuclide>;
class NReaction {
public:
  NReaction(const pugi::xml_node &reactionNode, NuclidePtr parent);
  std::string type_;
  std::string targetName_;
  std::shared_ptr<Nuclide> target_;
  std::shared_ptr<Nuclide> parent_;
  double branchingRatio_;
  double Q_;
};

#endif
