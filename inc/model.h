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

    explicit Model();

    void run();

    std::string chainpath() const;

    void set_chainpath(const std::string &chainpath);

    std::string inputpath() const;

    std::map<std::string, double> initcc_;
    std::string resultpath_;
    std::string solvertype_;
    std::vector<double> times_;

private:
    void readTimes(const YAML::Node &input);

    void readCc(const YAML::Node &input);

    void readSolverType(const YAML::Node &input);

    std::vector<double> linspace(const std::vector<std::string> &splat) const;

    std::vector<double> logspace(const std::vector<std::string> &splat) const;

    std::string chainpath_;

    Chain chain_;


    std::string inputpath_;


    const std::map<std::string, double> time_units = {
        {"s", 1.},
        {"h", 3600.},
        {"d", 3600. * 24.},
        {"y", 3600. * 24. * 365.25}
    };
};

#endif
