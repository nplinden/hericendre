#include "decay.h"
#include "nuclide.h"

Decay::Decay(const pugi::xml_node &decayNode, NuclidePtr parent) {
  parentName_ = parent->name_;
  parent_ = parent;

  type_ = decayNode.attribute("type").value();
  if (decayNode.attribute("target")) {
    targetName_ = decayNode.attribute("target").value();
    hasTarget_ = true;
    if (targetName_ == parentName_) {
      targetName_ = "";
      hasTarget_ = false;
    }
  } else {
    targetName_ = "";
    hasTarget_ = false;
  }

  std::string br = decayNode.attribute("branching_ratio").value();
  if (!br.empty())
    branchingRatio_ = stod(br);
  else
    branchingRatio_ = 1.;
};

Decay::Decay(std::string type, std::string targetName, double branchingRatio,
             NuclidePtr parent) {
  this->type_ = type;
  this->targetName_ = targetName;
  this->branchingRatio_ = branchingRatio;
  this->parentName_ = parent->name_;
  this->parent_ = parent;
}

const std::map<std::string, std::vector<std::string>> Decay::SECONDARIES = {
    {"alpha", {"He4"}},
    {"beta-,alpha", {"He4"}},
    {"ec/beta+,alpha", {"He4"}},
    {"ec/beta+,p", {"H1"}},
    {"ec/beta+,p,p", {"H1", "H1"}},
    {"ec/beta+,p,p,p", {"H1", "H1", "H1"}},
    {"p", {"H1"}},
    {"p,p", {"H1", "H1"}},
};
