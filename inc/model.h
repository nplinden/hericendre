#ifndef INPUT_HPP_INCLUDE
#define INPUT_HPP_INCLUDE
#include <map>
#include <string>
#include <vector>
#include "yaml-cpp/yaml.h"
#include "chain.h"

class Model {
public:
    // CONSTRUCTORS
    explicit Model(const std::string &inputpath);

    void run();

    std::string solvertype_;
    std::string result_path_;
    std::map<std::string, double> concentrations_;

    std::string chainpath() const;

    void set_chainpath(const std::string &chainpath);

    std::vector<double> times() const;

    void set_times(const std::vector<double> &times);

    std::string inputpath() const;

private:
    void readTimes(const YAML::Node &input);

    void readCc(const YAML::Node &input);

    void readSolverType(const YAML::Node &input);

    std::vector<double> linspace(const std::vector<std::string> &splat) const;

    std::vector<double> logspace(const std::vector<std::string> &splat) const;

    std::string chainpath_;
    Chain chain_;

    std::vector<double> times_;

    std::string inputpath_;

    const std::map<std::string, double> time_units = {
        {"s", 1.},
        {"h", 3600.},
        {"d", 3600. * 24.},
        {"y", 3600. * 24. * 365.25}
    };
};

#endif
