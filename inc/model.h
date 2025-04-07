#ifndef INPUT_HPP_INCLUDE
#define INPUT_HPP_INCLUDE
#include <map>
#include <string>
#include <vector>
#include "chain.h"
#include <toml++/toml.hpp>

enum SolverType
{
    DECAY,
    CRAM48
};

const std::map<std::string, double> TIME_UNITS = {
    {"s", 1.},
    {"h", 3600.},
    {"d", 3600. * 24.},
    {"y", 3600. * 24. * 365.25}};

class Model
{
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
    void readSettings(const toml::table &tbl);
    void readTime(const toml::table &tbl);
    std::vector<double> compute_time_function(const std::string &str);
    void readMaterial(const toml::table &tbl);

    std::vector<double> linspace(const std::vector<std::string> &splat) const;

    std::vector<double> logspace(const std::vector<std::string> &splat) const;

    std::string chainpath_;

    Chain chain_;

    std::string inputpath_;
    std::string name_;

    std::string time_unit_name_;
    double time_unit_magnitude_;
};

#endif
