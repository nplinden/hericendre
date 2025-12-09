#include "decay.h"
#include "nuclide.h"

Decay::Decay(const pugi::xml_node &decayNode, const NuclidePtr &parent) {
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

    if (const std::string br = decayNode.attribute("branching_ratio").value(); !br.empty())
        branchingRatio_ = stod(br);
    else
        branchingRatio_ = 1.;
};

Decay::Decay(const std::string &type, const std::string &targetName, const double &branchingRatio,
             const NuclidePtr &parent) {
    this->type_ = type;
    this->targetName_ = targetName;
    this->hasTarget_ = true;
    this->branchingRatio_ = branchingRatio;
    this->parentName_ = parent->name_;
    this->parent_ = parent;
}



Decay Decay::getSecondaries() const {
    if (this->hasSecondaries()) {
        auto secondary = SECONDARIES.at(this->type_);
        Decay secondaryDecay(
            this->type_,
            secondary.particle_,
            this->branchingRatio_ * secondary.multiplicity_,
            this->parent_);
        return secondaryDecay;
    }
    throw std::runtime_error("No secondaries for this decay type");
}

const std::map<std::string, SecondaryParticle> Decay::SECONDARIES = {
    {"alpha", SecondaryParticle("He4", 1)},
    {"beta-,alpha", SecondaryParticle("He4", 1)},
    {"ec/beta+,alpha", SecondaryParticle("He4", 1)},
    {"ec/beta+,p", SecondaryParticle("H1", 1)},
    {"ec/beta+,p,p", SecondaryParticle("H1", 2)},
    {"ec/beta+,p,p,p", SecondaryParticle("H1", 3)},
    {"p", SecondaryParticle("H1", 1)},
    {"p,p", SecondaryParticle("H1", 2)},
};
