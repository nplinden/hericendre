#ifndef DECAY_HPP_INCLUDED
#define DECAY_HPP_INCLUDED
#include <memory>
#include <pugixml.hpp>
#include <string>

class Nuclide;
using NuclidePtr = std::shared_ptr<Nuclide>;
class Decay {
public:
  // Constructor
  Decay(const pugi::xml_node &decayNode, NuclidePtr parent);

  // Member Variables
  std::string type_;
  std::string targetName_;
  std::string parentName_;
  std::shared_ptr<Nuclide> target_;
  std::shared_ptr<Nuclide> parent_;
  double branchingRatio_;
};

#endif
