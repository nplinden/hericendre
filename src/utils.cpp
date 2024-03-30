#include "utils.hpp"
#include <fmt/format.h>

std::vector<std::string> split(std::string str) {
  std::vector<std::string> vect;
  vect.push_back(std::string());
  for (const auto c : str) {
    if (c == ' ') {
      vect.push_back(std::string());
    } else {
      vect.back() += c;
    }
  }
  return vect;
}

std::string concatenate(std::vector<std::string> strs) {
  std::string tmp = "";
  for (const auto c : strs) {
    tmp += c;
    tmp += " ";
  }
  tmp.pop_back();
  return tmp;
}

std::string fmtDouble(double d) {
  std::string s = fmt::format("{}", d);
  if (d == floor(d))
    s += ".0";
  return s;
}
