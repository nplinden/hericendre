#include "reaction.h"
#include "nuclide.h"

Reaction::Reaction(const pugi::xml_node &reactionNode, const NuclidePtr &parent)
{
    parent_ = parent;
    type_ = reactionNode.attribute("type").value();

    if (reactionNode.attribute("Q"))
        Q_ = reactionNode.attribute("Q").as_double();

    if (reactionNode.attribute("target"))
    {
        targetName_ = reactionNode.attribute("target").value();
    }
    else
    {
        targetName_ = "";
    }

    branchingRatio_ = 1.0;
    if (reactionNode.attribute("branching_ratio"))
        branchingRatio_ = reactionNode.attribute("branching_ratio").as_double();
}