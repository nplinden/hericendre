#include "nfy.hpp"
#include "utils.hpp"
#include <iostream>
#include <fmt/core.h>
#include <fmt/os.h>
#include <fmt/ranges.h>

NeutronFissionYield::NeutronFissionYield(const pugi::xml_node& nfyNode){
    auto energies = split(std::string(nfyNode.child_value("energies"))) ;
    for (auto e: energies){
        fmt::print("{} ", e) ;
    }
    fmt::print("\n") ;
    
    for (pugi::xml_node fyNode = nfyNode.child("fission_yields"); fyNode; fyNode = fyNode.next_sibling("")){
        std::string energy = fyNode.attribute("energy").value() ;
        fy_[energy] = std::map<std::string, double>() ;
        std::vector<std::string> products = split(fyNode.child_value("products")) ;
        std::vector<std::string> data = split(fyNode.child_value("data")) ;
        for (int i = 0; i < products.size(); i++){
            fy_[energy][products[i]] = stod(data[i]) ;
        }
    }
}