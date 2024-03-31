#ifndef UTILS_HPP_INCLUDE
#define UTILS_HPP_INCLUDE
#include <string>
#include <vector>

/**
 * \brief Split a character string along spaces.
 *
 * \param str
 */
std::vector<std::string> split(std::string str);

/**
 * \brief Concatenate a vector of strings in a single string,
 * with space separators.
 *
 * \param A vector of strings
 */
std::string concatenate(std::vector<std::string> strs);

/**
 * \brief Format a double using fmt, if the double has an integer
 * value, add a .0 suffix.
 *
 *\param A double.
 */
std::string fmtDouble(double d);

#endif
