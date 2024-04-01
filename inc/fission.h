#ifndef NFY_HPP_INCLUDED
#define NFY_HPP_INCLUDED
#include <map>
#include <memory>
#include <nuclide.h>
#include <pugixml.hpp>
#include <string>
#include <vector>

using NuclidePtr = std::shared_ptr<Nuclide>;
class Fission {
public:
  /**
   * \brief Construct a new Fission object from a fission yield node and a
   * pointer to fissile parent.
   *
   * \param NFYNode
   * \param parent
   */
  Fission(const pugi::xml_node &NFYNode, NuclidePtr parent);

  /**
   * \brief Add a fission yield node to a parent node.
   *
   * \param nfyNode
   */
  void addNode(pugi::xml_node &nfyNode);

  /**
   * \brief The list of energies for which fission yields are defined.
   *
   */
  std::vector<std::string> energies_;

  /**
   * \brief The map containing fission yields values. First layer is the
   * incident neutron energy value as a string. The second layer is the
   * fission product name.
   *
   */
  std::map<std::string, std::map<std::string, double>> fy_;

  /**
   * \brief The name of the fissile nuclide.
   *
   */
  std::string parentName_;

  /**
   * \brief A pointer to the fissile nuclide.
   *
   */
  NuclidePtr parent_;
};

#endif
