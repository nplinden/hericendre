#include "utils.hpp"


std::vector<std::string> split(std::string str){
    std::vector<std::string> vect ;
    vect.push_back(std::string()) ;
    for (const auto c: str){
        if (c == ' '){
            vect.push_back(std::string()) ;
        } else {
            vect.back() += c ;
        }
    }
    return vect ;
}