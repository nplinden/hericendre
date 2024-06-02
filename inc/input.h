#ifndef INPUT_HPP_INCLUDE
#define INPUT_HPP_INCLUDE
#include <map>
#include <string>
#include <vector>
#include "yaml-cpp/yaml.h"
#include "chain.h"

class Input {
public:
  // CONSTRUCTORS
  explicit Input(const std::string &inputpath);

  void readTimes(const YAML::Node& input);
  void readCc(const YAML::Node& input);
  void readSolverType(const YAML::Node& input);
  void run();

  std::string inputpath_;
  std::string solvertype_ ;
  std::string chainpath_;
  Chain chain_ ;
  std::string result_path_ ;
  std::vector<double> times_;
  std::map<std::string, double> concentrations_;
  //
  // void addNode(pugi::xml_node &rootnode);
  std::vector<double> linspace(const std::vector<std::string>& splat) const;
  std::vector<double> logspace(const std::vector<std::string>& splat) const;
  std::vector<double> defaultline(const std::vector<std::string>& splat) const;
  const std::map<std::string, double> time_units = {{"s", 1.},
    {"h", 3600.},
    {"d", 3600.*24.},
    {"y", 3600.*24.*365.25}
  };
};

#endif
