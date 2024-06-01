#ifndef INPUT_HPP_INCLUDE
#define INPUT_HPP_INCLUDE
#include <fmt/core.h>
#include <map>
#include <string>
#include <vector>

class Input {
public:
  // CONSTRUCTORS
  explicit Input(const std::string &inputpath);

  void run() const;

  std::string inputpath_;
  std::string chainpath_;
  std::string cyclemode_;
  std::vector<double> times_;
  std::vector<double> powerlevel_;
  std::map<std::string, double> concentrations_;
  //
  // void addNode(pugi::xml_node &rootnode);
};

#endif
