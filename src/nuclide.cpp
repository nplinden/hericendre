#include "nuclide.hpp"
#include <algorithm>
#include <iostream>
#include <cmath>
#include <fmt/format.h>
#include "decay.hpp"
#include "nreaction.hpp"
#include "utils.hpp"

Nuclide::Nuclide(std::string name, double dconst){
    name_ = name ;
    dconst_ = dconst ;
} ;

Nuclide::Nuclide(const pugi::xml_node& nuclideNode){
    name_ = nuclideNode.attribute("name").value() ;
    std::string halflife = nuclideNode.attribute("half_life").value() ;
    dconst_ = !halflife.empty() ? std::log(2) / stod(halflife) : 0. ;

    std::tuple<int, int, int> zam_tuple = getZam(name_) ;
    z_ = std::get<0>(zam_tuple) ;
    a_ = std::get<1>(zam_tuple) ;
    m_ = std::get<2>(zam_tuple) ;
    zam_ = 10000 * z_ + 10 * a_ + m_ ;

    std::string nreac_str = nuclideNode.attribute("reactions").value() ;
    nreac_ = !nreac_str.empty() ? stoi(nreac_str) : 0 ;

    std::string decay_str = nuclideNode.attribute("decay_modes").value() ;
    ndecay_ = !decay_str.empty() ? stoi(decay_str) : 0 ;

    std::string denergy_str = nuclideNode.attribute("decay_energy").value() ;
    denergy_ = !denergy_str.empty() ? stod(denergy_str) : 0. ;

    for (pugi::xml_node sourceNode = nuclideNode.child("source"); sourceNode; sourceNode = sourceNode.next_sibling("source")){
        sources_.push_back(std::make_shared<Source>(Source(sourceNode))) ;
    }

    // auto nfyNode = nuclideNode.child("neutron_fission_yields") ;
    // if (nfyNode){
    //     std::string parent = nfyNode.attribute("parent").value() ;
    //     if (!parent.empty()) {
    //         this->nfyParent_ = parent ;
    //     } else if (nfyNode.child("energies")){
    //         neutronFissionYields_ = std::make_shared<NeutronFissionYield>(nfyNode) ;
    //     }
    // }
} ;

const std::map<std::string, int> Nuclide::ELEMENTS = {
    {"n", 0}, {"H", 1}, {"He", 2}, {"Li", 3}, {"Be", 4}, {"B", 5}, {"C", 6}, {"N", 7}, {"O", 8}, {"F", 9}, {"Ne", 10}, 
    {"Na", 11}, {"Mg", 12}, {"Al", 13}, {"Si", 14}, {"P", 15}, {"S", 16}, {"Cl", 17}, {"Ar", 18}, {"K", 19}, {"Ca", 20}, 
    {"Sc", 21}, {"Ti", 22}, {"V", 23}, {"Cr", 24}, {"Mn", 25}, {"Fe", 26}, {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30}, 
    {"Ga", 31}, {"Ge", 32}, {"As", 33}, {"Se", 34}, {"Br", 35}, {"Kr", 36}, {"Rb", 37}, {"Sr", 38}, {"Y", 39}, {"Zr", 40}, 
    {"Nb", 41}, {"Mo", 42}, {"Tc", 43}, {"Ru", 44}, {"Rh", 45}, {"Pd", 46}, {"Ag", 47}, {"Cd", 48}, {"In", 49}, {"Sn", 50}, 
    {"Sb", 51}, {"Te", 52}, {"I", 53}, {"Xe", 54}, {"Cs", 55}, {"Ba", 56}, {"La", 57}, {"Ce", 58}, {"Pr", 59}, {"Nd", 60}, 
    {"Pm", 61}, {"Sm", 62}, {"Eu", 63}, {"Gd", 64}, {"Tb", 65}, {"Dy", 66}, {"Ho", 67}, {"Er", 68}, {"Tm", 69}, {"Yb", 70}, 
    {"Lu", 71}, {"Hf", 72}, {"Ta", 73}, {"W", 74}, {"Re", 75}, {"Os", 76}, {"Ir", 77}, {"Pt", 78}, {"Au", 79}, {"Hg", 80}, 
    {"Tl", 81}, {"Pb", 82}, {"Bi", 83}, {"Po", 84}, {"At", 85}, {"Rn", 86}, {"Fr", 87}, {"Ra", 88}, {"Ac", 89}, {"Th", 90}, 
    {"Pa", 91}, {"U", 92}, {"Np", 93}, {"Pu", 94}, {"Am", 95}, {"Cm", 96}, {"Bk", 97}, {"Cf", 98}, {"Es", 99}, {"Fm", 100}, 
    {"Md", 101}, {"No", 102}, {"Lr", 103}, {"Rf", 104}, {"Db", 105}, {"Sg", 106}, {"Bh", 107}, {"Hs", 108}, {"Mt", 109}, 
    {"Ds", 110}, {"Rg", 111}
    } ;

std::tuple<int, int, int> Nuclide::getZam(std::string name){
    std::string element ;
    std::string mass_str ;
    std::string meta_str ;
    bool elem_flag = true ;
    bool mass_flag = false ;
    for (char& c : name){
        if (!isdigit(c) & elem_flag){
            element += c ;
        } else if (isdigit(c) & elem_flag){
            elem_flag = false ;
            mass_flag = true ;
            mass_str += c ; 
        } else if (isdigit(c) & mass_flag){
            mass_str += c ; 
        } else {
            mass_flag = false ;
            meta_str += c ;
        }
    }
    std::string meta_number_str ;
    for (char& c : meta_str){
        if (isdigit(c)) meta_number_str += c ;
    }
    int z = ELEMENTS.at(element) ;
    int a = stoi(mass_str) ;
    int m = meta_number_str.empty() ? 0 : stoi(meta_number_str) ;
    return {z, a, m} ;
}

std::string Nuclide::str() {
    return name_ ;
}

bool Nuclide::operator<(const Nuclide& other) const {
    return zam_ < other.zam_ ;
}

void Nuclide::addNode(pugi::xml_node& rootnode){
    auto nucNode = rootnode.append_child("nuclide") ;
    nucNode.append_attribute("name") = this->name_.c_str() ;
    if (this->dconst_ != 0.)
        nucNode.append_attribute("half_life") = fmt::format("{}", std::log(2) / this->dconst_).c_str() ;
    if (this->ndecay_ != 0.)
        nucNode.append_attribute("decay_modes") = this->ndecay_ ;
    if (this->denergy_ != 0.)
        nucNode.append_attribute("decay_energy") = fmtDouble(this->denergy_).c_str() ;
    nucNode.append_attribute("reactions") = this->nreac_ ;

    if (!this->decays_.empty()){
        for (auto decay : this->decays_){
            auto decNode = nucNode.append_child("decay") ;
            if (!decay->type_.empty())
                decNode.append_attribute("type") = decay->type_.c_str();
            if (!decay->targetName_.empty())
                decNode.append_attribute("target") = decay->targetName_.c_str();
            decNode.append_attribute("branching_ratio") = fmtDouble(decay->branchingRatio_).c_str() ;
        }
    }

    if (!this->sources_.empty()){
        for (auto source : this->sources_){
            auto sourceNode = nucNode.append_child("source") ;
            sourceNode.append_attribute("type") = source->type_.c_str() ;
            if (!source->interpolation_.empty())
                sourceNode.append_attribute("interpolation") = source->interpolation_.c_str() ;
            sourceNode.append_attribute("particle") = source->particle_.c_str() ;
            if (source->type_ != "mixture"){
                std::string parameters = "" ;
                for (auto energy : source->energy_){
                    parameters += fmtDouble(energy) ;
                    parameters += " " ;
                }
                for (auto intensity : source->intensity_){
                    parameters += fmt::format("{} ", intensity) ;
                }
                parameters.pop_back() ;
                auto paramNode = sourceNode.append_child("parameters") ;
                paramNode.text() = parameters.c_str() ;
            } else {
                for (int i = 0; i < source->pairs_.size(); i++){
                    auto p = source->pairs_[i] ;
                    auto proba = source->pair_probabilities_[i] ;
                    auto pairNode = sourceNode.append_child("pair") ;
                    pairNode.append_attribute("probability") = proba ;
                    auto distNode = pairNode.append_child("dist") ;
                    distNode.append_attribute("type") = p.type_.c_str() ;
                    if (!p.interpolation_.empty())
                        distNode.append_attribute("interpolation") = p.interpolation_.c_str() ;

                    std::string parameters = "" ;
                    for (auto energy : p.energy_)
                        parameters += fmt::format("{} ", energy) ;
                    for (auto intensity : p.intensity_)
                        parameters += fmt::format("{} ", intensity) ;
                    parameters.pop_back() ;
                    auto paramNode = distNode.append_child("parameters") ;
                    paramNode.text() = parameters.c_str() ;
                }
            }
        }
    }

    if (!this->reactions_.empty()){
        for (auto reaction : this->reactions_){
            auto reacNode = nucNode.append_child("reaction") ;
            reacNode.append_attribute("type") = reaction->type_.c_str() ;
            reacNode.append_attribute("Q") = fmtDouble(reaction->Q_).c_str() ;
            if (!reaction->targetName_.empty())
                reacNode.append_attribute("target") = reaction->targetName_.c_str() ;
            if (reaction->branchingRatio_ != 1.)
                reacNode.append_attribute("branching_ratio") = reaction->branchingRatio_ ;
        }
    }

    // if (this->neutronFissionYields_){
    //     // fmt::print("{}", this->name_) ;
    //     auto nfyNode = nucNode.append_child("neutron_fission_yields") ;
    //     if (!this->nfyParent_.empty()){
    //         nfyNode.append_attribute("parent") = this->nfyParent_.c_str() ;
    //     } else {
    //         this->neutronFissionYields_->addNode(nfyNode) ;
    //     }
    // }
}