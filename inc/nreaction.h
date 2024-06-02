#ifndef NREACTION_HPP_INCLUDED
#define NREACTION_HPP_INCLUDED
#include <memory>
#include <pugixml.hpp>
#include <string>

class Nuclide;
using NuclidePtr = std::shared_ptr<Nuclide>;

class NReaction {
public:
 /**
  * \brief An object describing a neutronic reaction. In Hericendre a reaction
  * is defined as
  *
  *              (parent) -(type)[branching ratio]-> (target).
  *
  * If a (parent, type) pair has k branching ratios and targets, k NReaction
  * object will be needed.
  *
  *
  * \param
  */
 NReaction(const pugi::xml_node &reactionNode, NuclidePtr parent);

 /**
  * \brief The type of reaction. Possible types include (n,gamma), (n,a),
  * (n,p), (n,2n), etc.
  */
 std::string type_;
 std::string targetName_;

 /**
  * \brief A pointer to the reaction target object.
  */
 std::shared_ptr<Nuclide> target_;

 /**
  * \brief A pointer to the reaction parent object.
  */
 std::shared_ptr<Nuclide> parent_;

 /**
  * \brief The reaction's branching ratio.
  */
 double branchingRatio_;

 /**
  * \brief The reaction's Q value.
  */
 double Q_;
};

#endif
