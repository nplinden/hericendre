#ifndef MICROXS_HPP_INCLUDED
#define MICROXS_HPP_INCLUDED
#include <string>
#include <vector>

class MicroXS
{
public:
    MicroXS();
    MicroXS(const std::string &xspath);
    double getXS(const size_t nuclideIndex,
                 const size_t reactionIndex) const
    {
        return data_.at(nuclideIndex).at(reactionIndex);
    };

    double getXS(const std::string &nuclideName,
                 const std::string &reactionType) const;

private:
    std::vector<std::vector<double>> data_;
    std::vector<std::string> nuclides_;
    std::vector<std::string> reactions_;
};

#endif