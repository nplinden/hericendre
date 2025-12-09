#ifndef SECONDARIES_HPP_INCLUDED
#define SECONDARIES_HPP_INCLUDED
#include <string>

class SecondaryParticle {
    public:
    SecondaryParticle(const std::string &particle, int multiplicity): particle_(particle), multiplicity_(multiplicity) {}
    std::string particle_;
    int multiplicity_;
};

#endif // SECONDARIES_HPP_INCLUDED