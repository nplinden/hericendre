#ifndef UTILS_HPP_INCLUDE
#define UTILS_HPP_INCLUDE
#include <string>
#include <vector>
#include <map>

/**
 * \brief Split a character string along spaces.
 *
 * \param str
 */
std::vector<std::string> split(const std::string &str);

/**
 * \brief Trim whitespaces from a string.
 *
 * \param str
 */
std::string trim(const std::string &str);

/**
 * \brief Concatenate a vector of strings in a single string,
 * with space separators.
 *
 * \param strs: A vector of strings
 */
std::string concatenate(const std::vector<std::string> &strs);

/**
 * \brief Format a double using fmt, if the double has an integer
 * value, add a .0 suffix.
 *
 *\param d: A double.
 */
std::string fmtDouble(double d);

template <typename K, typename V>
V atWithDefault(const std::map<K, V> &m, const K &k, const V &deflt)
{
    if (m.find(k) != m.end())
    {
        return m.at(k);
    }
    return deflt;
    // if (m.contains(k)) return m.at(k) ;
    // return deflt;
};

#endif
