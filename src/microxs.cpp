#include <microxs.h>
#include <fmt/core.h>
#include <highfive/H5Easy.hpp>
#include <utils.h>

using Vector3D = std::vector<std::vector<std::vector<double>>>;
using MicroXSPtr = std::shared_ptr<MicroXS>;

MicroXS::MicroXS() = default;

MicroXS::MicroXS(const std::string &xspath)
{
    H5Easy::File file(xspath, H5Easy::File::ReadOnly);

    auto data = H5Easy::load<Vector3D>(file, "/data");
    for (const auto &nuclide_dim : data)
    {
        data_.push_back(std::vector<double>{});
        for (const auto &reaction_dim : nuclide_dim)
        {
            if (reaction_dim.size() != 1)
            {
                throw std::runtime_error("MicroXS currently only supports 1 group cross-section data.");
            }
            data_.back().push_back(reaction_dim[0]);
        }
    }
    for (auto &nuclide : H5Easy::load<std::vector<std::string>>(file, "/nuclides"))
    {
        nuclide.erase(nuclide.find_last_not_of('\0') + 1);
        nuclides_.push_back(nuclide);
    }

    for (auto &reaction : H5Easy::load<std::vector<std::string>>(file, "/reactions"))
    {
        reaction.erase(reaction.find_last_not_of('\0') + 1);
        reactions_.push_back(trim(reaction));
    }
}

double MicroXS::getXS(const std::string &nuclideName,
                      const std::string &reactionType) const
{
    auto nuclideIt = std::find(nuclides_.begin(), nuclides_.end(), nuclideName);
    if (nuclideIt == nuclides_.end())
    {
        throw std::runtime_error(fmt::format("Nuclide {} not found in MicroXS data.", nuclideName));
    }
    size_t nuclideIndex = std::distance(nuclides_.begin(), nuclideIt);

    auto reactionIt = std::find(reactions_.begin(), reactions_.end(), reactionType);
    if (reactionIt == reactions_.end())
    {
        throw std::runtime_error(fmt::format("Reaction type {} not found in MicroXS data.", reactionType));
    }
    size_t reactionIndex = std::distance(reactions_.begin(), reactionIt);

    return getXS(nuclideIndex, reactionIndex);
}