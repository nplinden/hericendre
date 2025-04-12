#include "utils.h"
#include <fmt/format.h>

std::vector<std::string> split(const std::string &str)
{
    std::vector<std::string> vect;
    vect.emplace_back();
    for (const auto c : str)
    {
        if (c == ' ')
        {
            vect.emplace_back();
        }
        else
        {
            vect.back() += c;
        }
    }
    return vect;
}

std::string trim(const std::string &str)
{
    size_t start = str.find_first_not_of(" \t\n\r\f\v");

    if (start == std::string::npos)
    {
        return "";
    }

    size_t end = str.find_last_not_of(" \t\n\r\f\v");
    return str.substr(start, end - start + 1);
}

std::string concatenate(const std::vector<std::string> &strs)
{
    std::string tmp;
    for (const auto &c : strs)
    {
        tmp += c;
        tmp += " ";
    }
    tmp.pop_back();
    return tmp;
}

std::string fmtDouble(double d)
{
    std::string s = fmt::format("{}", d);
    if (d == floor(d))
        s += ".0";
    return s;
}

// template <typename K, typename V>
// V atWithDefault(const std::map<K, V>& m, const K& k, const V& deflt){
//   if (m.contains(k)) return m.at(k) ;
//   return deflt;
// }
