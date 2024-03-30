#ifndef UTILS_HPP_INCLUDE
#define UTILS_HPP_INCLUDE
#include <string>
#include <vector>

std::vector<std::string> split(std::string str);
std::string concatenate(std::vector<std::string> strs);
std::string fmtDouble(double d);

#endif