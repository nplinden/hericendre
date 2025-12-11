#include "nuclide.h"
#include "decay.h"
#include "utils.h"
#include <cmath>
#include <fmt/format.h>

Nuclide::Nuclide(const std::string &name, const double dconst)
{
    name_ = name;
    dconst_ = dconst;
};

Nuclide::Nuclide(const pugi::xml_node &nuclideNode)
{
    name_ = nuclideNode.attribute("name").value();
    const std::string halflife = nuclideNode.attribute("half_life").value();
    dconst_ = !halflife.empty() ? std::log(2) / stod(halflife) : 0.;

    const std::tuple<int, int, int> zam_tuple = getZam(name_);
    z_ = std::get<0>(zam_tuple);
    a_ = std::get<1>(zam_tuple);
    m_ = std::get<2>(zam_tuple);
    zam_ = 10000 * z_ + 10 * a_ + m_;

    const std::string nreac_str = nuclideNode.attribute("reactions").value();
    nreac_ = !nreac_str.empty() ? stoi(nreac_str) : 0;

    const std::string decay_str = nuclideNode.attribute("decay_modes").value();
    ndecay_ = !decay_str.empty() ? stoi(decay_str) : 0;

    const std::string denergy_str = nuclideNode.attribute("decay_energy").value();
    denergy_ = !denergy_str.empty() ? stod(denergy_str) : 0.;
};

const std::map<std::string, int> Nuclide::ELEMENTS = {
    {"n", 0}, {"H", 1}, {"He", 2}, {"Li", 3}, {"Be", 4}, {"B", 5}, {"C", 6}, {"N", 7}, {"O", 8}, {"F", 9}, {"Ne", 10}, {"Na", 11}, {"Mg", 12}, {"Al", 13}, {"Si", 14}, {"P", 15}, {"S", 16}, {"Cl", 17}, {"Ar", 18}, {"K", 19}, {"Ca", 20}, {"Sc", 21}, {"Ti", 22}, {"V", 23}, {"Cr", 24}, {"Mn", 25}, {"Fe", 26}, {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30}, {"Ga", 31}, {"Ge", 32}, {"As", 33}, {"Se", 34}, {"Br", 35}, {"Kr", 36}, {"Rb", 37}, {"Sr", 38}, {"Y", 39}, {"Zr", 40}, {"Nb", 41}, {"Mo", 42}, {"Tc", 43}, {"Ru", 44}, {"Rh", 45}, {"Pd", 46}, {"Ag", 47}, {"Cd", 48}, {"In", 49}, {"Sn", 50}, {"Sb", 51}, {"Te", 52}, {"I", 53}, {"Xe", 54}, {"Cs", 55}, {"Ba", 56}, {"La", 57}, {"Ce", 58}, {"Pr", 59}, {"Nd", 60}, {"Pm", 61}, {"Sm", 62}, {"Eu", 63}, {"Gd", 64}, {"Tb", 65}, {"Dy", 66}, {"Ho", 67}, {"Er", 68}, {"Tm", 69}, {"Yb", 70}, {"Lu", 71}, {"Hf", 72}, {"Ta", 73}, {"W", 74}, {"Re", 75}, {"Os", 76}, {"Ir", 77}, {"Pt", 78}, {"Au", 79}, {"Hg", 80}, {"Tl", 81}, {"Pb", 82}, {"Bi", 83}, {"Po", 84}, {"At", 85}, {"Rn", 86}, {"Fr", 87}, {"Ra", 88}, {"Ac", 89}, {"Th", 90}, {"Pa", 91}, {"U", 92}, {"Np", 93}, {"Pu", 94}, {"Am", 95}, {"Cm", 96}, {"Bk", 97}, {"Cf", 98}, {"Es", 99}, {"Fm", 100}, {"Md", 101}, {"No", 102}, {"Lr", 103}, {"Rf", 104}, {"Db", 105}, {"Sg", 106}, {"Bh", 107}, {"Hs", 108}, {"Mt", 109}, {"Ds", 110}, {"Rg", 111}};

std::tuple<int, int, int> Nuclide::getZam(const std::string &name)
{
    std::string element;
    std::string mass_str;
    std::string meta_str;
    bool elem_flag = true;
    bool mass_flag = false;
    for (const char &c : name)
    {
        if (!isdigit(c) & elem_flag)
        {
            element += c;
        }
        else if (isdigit(c) & elem_flag)
        {
            elem_flag = false;
            mass_flag = true;
            mass_str += c;
        }
        else if (isdigit(c) & mass_flag)
        {
            mass_str += c;
        }
        else
        {
            mass_flag = false;
            meta_str += c;
        }
    }
    std::string meta_number_str;
    for (const char &c : meta_str)
    {
        if (isdigit(c))
            meta_number_str += c;
    }
    int z = ELEMENTS.at(element);
    int a = stoi(mass_str);
    int m = meta_number_str.empty() ? 0 : stoi(meta_number_str);
    return {z, a, m};
}

std::string Nuclide::str() const { return name_; }

bool Nuclide::operator<(const Nuclide &other) const
{
    return zam_ < other.zam_;
}

void Nuclide::addNode(pugi::xml_node &rootnode)
{
    auto nucNode = rootnode.append_child("nuclide");
    nucNode.append_attribute("name") = this->name_.c_str();
    if (this->dconst_ != 0.)
        nucNode.append_attribute("half_life") =
            fmt::format("{}", std::log(2) / this->dconst_).c_str();
    if (this->ndecay_ != 0.)
        nucNode.append_attribute("decay_modes") = this->ndecay_;
    if (this->denergy_ != 0.)
        nucNode.append_attribute("decay_energy") =
            fmtDouble(this->denergy_).c_str();
    nucNode.append_attribute("reactions") = this->nreac_;

    if (!this->decays_.empty())
    {
        for (const auto &decay : this->decays_)
        {
            auto decNode = nucNode.append_child("decay");
            if (!decay->type_.empty())
                decNode.append_attribute("type") = decay->type_.c_str();
            if (!decay->targetName_.empty())
                decNode.append_attribute("target") = decay->targetName_.c_str();
            decNode.append_attribute("branching_ratio") =
                fmtDouble(decay->branchingRatio_).c_str();
        }
    }
}

bool Nuclide::isStable()
{
    return dconst_ == 0.;
}