#include "nreaction.hpp"

NReaction::NReaction(const pugi::xml_node &reactionNode, NuclidePtr parent) {
  type_ = reactionNode.attribute("type").value();
  targetName_ = reactionNode.attribute("target").value();

  std::string br = reactionNode.attribute("branching_ratio").value();
  if (!br.empty())
    branchingRatio_ = stod(br);
  else
    branchingRatio_ = 1.;

  std::string q = reactionNode.attribute("Q").value();
  if (!q.empty())
    Q_ = stod(q);
  else
    Q_ = 0.;

  parent_ = parent;
}
