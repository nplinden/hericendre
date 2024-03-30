#ifndef NFY_HPP_INCLUDED
#define NFY_HPP_INCLUDED
#include <map>
#include <memory>
#include <nuclide.hpp>
#include <pugixml.hpp>
#include <string>
#include <vector>

using NuclidePtr = std::shared_ptr<Nuclide>;
class Fission {
public:
  /**
   * @brief Construct a new Fission object from a fission yield node and a
   * pointer to fissile parent.
   *
   * @param NFYNode
   * @param parent
   */
  Fission(const pugi::xml_node &NFYNode, NuclidePtr parent);

  /**
   * @brief Add a fission yield node to a parent node.
   *
   * @param nfyNode
   */
  void addNode(pugi::xml_node &nfyNode);

  /**
   * @brief
   *
   */
  std::vector<std::string> energies_;
  std::map<std::string, std::map<std::string, double>> fy_;
  std::string parentName_;
  NuclidePtr parent_;
};

#endif
