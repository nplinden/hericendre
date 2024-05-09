#ifndef DECAY_HPP_INCLUDED
#define DECAY_HPP_INCLUDED
#include <map>
#include <memory>
#include <pugixml.hpp>
#include <string>
#include <vector>

class Nuclide;
using NuclidePtr = std::shared_ptr<Nuclide>;
class Decay {
public:
  // CONSTRUCTOR
  /**
   * \brief Creates a Decay object using an decay-tagged xml node and a pointer
   * to the decaying nuclide.
   *
   */
  Decay(const pugi::xml_node &decayNode, NuclidePtr parent);

  Decay(std::string type, std::string targetName_, double branchingRatio_,
        NuclidePtr parent);

  // MEMBER VARIABLES
  /**
   * \brief Type of decay. Can by sf, alpha, IT, beta-, ec/beta+, p, n, etc.
   */
  std::string type_;

  /**
   * \brief Name of the resulting nuclide.
   */
  std::string targetName_;

  /**
   * \brief Name of the decaying nuclide.
   */
  std::string parentName_;

  /**
   * \brief A pointer to the resulting nuclide.
   */
  std::shared_ptr<Nuclide> target_;

  /**
   * \brief A pointer to the parent nuclide.
   */
  std::shared_ptr<Nuclide> parent_;

  /**
   * \brief A branching value for the decay reaction.
   */
  double branchingRatio_;

  static const std::map<std::string, std::vector<std::string>> SECONDARIES;
};

#endif
