#ifndef UTILS_HPP_INCLUDE
#define UTILS_HPP_INCLUDE
#include <map>
#include <string>
#include <vector>

/**
 * @brief Splits a string into a vector of substrings based on space characters
 *
 * This function splits the input string at each space character and returns
 * a vector containing the separated substrings. Empty strings are included
 * if there are consecutive spaces.
 *
 * @param str The input string to split
 * @return std::vector<std::string> A vector containing the split substrings
 */
std::vector<std::string> split(const std::string &str);

/**
 * @brief Removes leading and trailing whitespace from a string
 *
 * This function removes all leading and trailing whitespace characters
 * including spaces, tabs, newlines, carriage returns, form feeds, and vertical
 * tabs.
 *
 * @param str The input string to trim
 * @return std::string A new string with leading and trailing whitespace
 * removed. Returns an empty string if the input contains only whitespace.
 */
std::string trim(const std::string &str);

/**
 * @brief Concatenates a vector of strings into a single space-separated string
 *
 * This function joins all strings in the input vector with a single space
 * character between each element.
 *
 * @param strs The vector of strings to concatenate
 * @return std::string A single string with all input strings joined by spaces
 */
std::string concatenate(const std::vector<std::string> &strs);

/**
 * @brief Formats a double value as a string, ensuring decimals are shown
 *
 * This function converts a double to a string using fmt::format. If the value
 * is a whole number (equal to its floor), it appends ".0" to ensure the decimal
 * point is visible in the output.
 *
 * @param d The double value to format
 * @return std::string The formatted string representation of the double
 */
std::string fmtDouble(double d);

/**
 * @brief Retrieves a value from a map with a default fallback
 *
 * This template function safely accesses a value in a map by key. If the key
 * exists, it returns the associated value. If the key is not found, it returns
 * the provided default value instead of throwing an exception.
 *
 * @tparam K The key type
 * @tparam V The value type
 * @param m The map to search in
 * @param k The key to look up
 * @param deflt The default value to return if the key is not found
 * @return V The value associated with the key, or the default value if key not
 * found
 */
template <typename K, typename V>
V atWithDefault(const std::map<K, V> &m, const K &k, const V &deflt) {
  if (m.find(k) != m.end()) {
    return m.at(k);
  }
  return deflt;
};

#endif
