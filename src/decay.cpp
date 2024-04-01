#include "decay.h"
#include "nuclide.h"

Decay::Decay(const pugi::xml_node &decayNode, NuclidePtr parent) {
  type_ = decayNode.attribute("type").value();
  targetName_ = decayNode.attribute("target").value();

  std::string br = decayNode.attribute("branching_ratio").value();
  if (!br.empty())
    branchingRatio_ = stod(br);
  else
    branchingRatio_ = 1.;

  parentName_ = parent->name_;
  parent_ = parent;
};
