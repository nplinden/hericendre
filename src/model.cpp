#include "model.h"
#include "chain.h"
#include "results.h"
#include "cramsolver.h"
#include "utils.h"
#include <cmath>
#include "yaml-cpp/yaml.h"
#include <fmt/ranges.h>
#include <fmt/core.h>
#include <algorithm>
#include <decaysolver.h>
#include <highfive/H5Easy.hpp>
#include <toml++/toml.hpp>
#include <optional>

Model::Model(const std::string &inputpath)
{
    toml::table tbl;

    try
    {
        tbl = toml::parse_file(inputpath);
    }
    catch (const toml::parse_error &err)
    {
        fmt::print("[ERROR] Parsing input file failed\n");
        throw std::runtime_error(err);
    }

    std::optional<std::string> name = tbl["name"].value<std::string>();
    this->name_ = name ? *name : "";

    readSettings(tbl);
    readTimes(tbl);

    /*
    readCc(input);
    */
}

Model::Model() = default;

void Model::readSettings(const toml::table &tbl)
{
    std::optional<std::string> chainpath = tbl["Settings"]["chain"].value<std::string>();
    if (chainpath)
    {
        chainpath_ = *chainpath;
        this->chain_ = Chain(chainpath_.c_str());
    }
    else
    {
        fmt::print("[ERROR] No depletion chain file was specified\n");
        throw std::runtime_error("No depletion chain file was specified");
    }

    std::optional<std::string> resultpath = tbl["Settings"]["results"].value<std::string>();
    if (resultpath)
    {
        this->resultpath_ = *resultpath;
    }
    else
    {
        fmt::print("[ERROR] No result file was specified\n");
        throw std::runtime_error("No result file was specified");
    }

    std::optional<std::string> solver = tbl["Settings"]["solver"].value<std::string>();
    std::set<std::string> allowed_solvers = {"CRAM48", "Decay"};
    if (!solver)
    {
        fmt::print("No solver type specified, using CRAM48\n");
        this->solvertype_ = "CRAM48";
    }
    else if (allowed_solvers.find(*solver) != allowed_solvers.end())
    {
        this->solvertype_ = *solver;
    }
    else
    {
        fmt::print("[ERROR] Invalid solver type \"{}\", allowed solver: {}\n", *solver, fmt::join(allowed_solvers, ", "));
        throw std::runtime_error("Invalid solver type");
    }
}

void Model::readTimes(const toml::table &tbl)
{
    // TIMESTAMPS
    auto timestamps = tbl["Time"]["timestamps"].as_array();
    if (!timestamps)
    {
        fmt::print("[ERROR] No timestamps provided\n");
        throw std::runtime_error("No timestamps provided");
    }
    else
    {
        for (const auto &item : *timestamps)
        {
            if (item.is<double>())
            {
                this->times_.push_back(*item.value<double>());
            }
            else if (item.is<std::string>())
            {
                auto expended = this->compute_time_function(*item.value<std::string>());
                for (const auto &t : expended)
                {
                    this->times_.push_back(t);
                }
            }
            else if (item.is<int64_t>())
            {
                this->times_.push_back(*item.value<double>());
            }
            else
            {
                fmt::print("[ERROR] Invalid type in timestamps definition\n");
                throw std::runtime_error("Invalid type in timestamps definition");
            }
        }
    }

    auto units = tbl["Time"]["unit"];
    if (!units)
    {
        this->time_unit_name_ = "s";
        this->time_unit_magnitude_ = 1.;
    }
    else if (units.is<std::string>())
    {
        auto unit_name = *units.value<std::string>();
        if (TIME_UNITS.find(unit_name) == TIME_UNITS.end())
        {
            fmt::print("[ERROR] Invalid unit name\n");
            throw std::runtime_error("Invalid unit name");
        }
        else
        {
            this->time_unit_name_ = unit_name;
            this->time_unit_magnitude_ = TIME_UNITS.at(unit_name);
        }
    }
    else if (units.is<toml::table>())
    {
        auto unit_table = *units.as_table();
        auto name = unit_table["name"];
        if (!name || !name.is<std::string>())
        {
            fmt::print("[ERROR] Invalid time unit name.\n");
            throw std::runtime_error("[ERROR] Invalid time unit name.");
        }
        else
        {
            this->time_unit_name_ = *name.value<std::string>();
        }
        auto magnitude = unit_table["magnitude"];
        if (!magnitude)
        {
            fmt::print("[ERROR] Invalid time unit magnitude.\n");
            throw std::runtime_error("[ERROR] Invalid time unit magnitude.");
        }
        else if (magnitude.is<double>() || magnitude.is<int64_t>())
        {
            this->time_unit_magnitude_ = *magnitude.value<double>();
        }
        else
        {
            fmt::print("[ERROR] Invalid time unit magnitude.\n");
            throw std::runtime_error("[ERROR] Invalid time unit magnitude.");
        }
    }
    else
    {
        fmt::print("[ERROR] Time unit must be either a string or a table.\n");
        throw std::runtime_error("Time unit must be either a string or a table.");
    }

    for (auto &t : this->times_)
    {
        t *= this->time_unit_magnitude_;
    }

    if (!std::is_sorted(this->times_.begin(), this->times_.end()))
    {
        fmt::print("[ERROR] Time vector is not sorted.\n");
        throw std::runtime_error("Time vector is not sorted.");
    }

    for (const double &t : this->times_)
    {
        fmt::print("{}\n", t);
    }
}

void Model::readCc(const YAML::Node &input)
{
    std::string ccmode = "Explicit";
    if (const auto mode = input["ConcentrationMode"])
        ccmode = mode.as<std::string>();
    if (ccmode == "Explicit")
    {
        const std::vector<std::string> cc = split(input["Concentrations"].as<std::string>());
        for (size_t i = 0; i < cc.size(); i += 2)
            initcc_[cc[i]] = stod(cc[i + 1]);
    }
    else if (ccmode == "Uniform")
    {
        const double val = stod(input["Concentrations"].as<std::string>());
        for (const auto &nuclide : chain_.nuclides_)
            initcc_[nuclide->name_] = val;
    }
    else
    {
        fmt::print("[ERROR] Invalid ConcentrationMode '{:s}'\n", ccmode);
        throw std::runtime_error(fmt::format("Invalid ConcentrationMode '{}'", ccmode));
    }
}

std::vector<double> Model::compute_time_function(const std::string &str)
{
    std::vector<std::string> splat = split(str);
    if (splat.size() < 4)
    {
        fmt::print("[ERROR] Invalid time function in timestamps definition");
        throw std::runtime_error("Invalid time function in timestamps definition");
    }
    else
    {
        auto kind = splat[0];
        if (kind == "linspace")
        {
            return this->linspace(splat);
        }
        else if (kind == "logspace")
        {
            return this->logspace(splat);
        }
        else
        {
            fmt::print("[ERROR] Invalid time function in timestamps definition");
            throw std::runtime_error("Invalid time function in timestamps definition");
        }
    }
}

std::vector<double> Model::linspace(const std::vector<std::string> &splat) const
{
    const double start = stod(splat[1]);
    const double stop = stod(splat[2]);
    const int nstep = stoi(splat[3]);
    const std::string unit = (splat.size() == 5) ? splat[4] : "s";

    const double step = (stop - start) / nstep;
    std::vector<double> values;
    for (int istep = 0; istep < nstep; istep++)
    {
        values.push_back((start + istep * step));
    }
    return values;
}

std::vector<double> Model::logspace(const std::vector<std::string> &splat) const
{
    const double start = std::log10(stod(splat[1]));
    const double stop = std::log10(stod(splat[2]));
    const int nstep = stoi(splat[3]);
    const std::string unit = (splat.size() == 5) ? splat[4] : "s";

    const double step = (stop - start) / nstep;
    std::vector<double> values;
    for (int istep = 0; istep < nstep; istep++)
    {
        values.push_back((pow(10, start + istep * step)));
    }
    return values;
}

void Model::run()
{
    if (solvertype_ == "Decay")
    {
        DecaySolver solver(chain_);
        solver.run(initcc_, times_);
        auto results = solver.results_;
        results.to_csv(resultpath_);
        H5Easy::File file("results.h5", H5Easy::File::Overwrite);
        results.to_hdf5(file);
        solver.to_hdf5((file));
    }
    else if (solvertype_ == "CRAM48")
    {
        CRAMSolver solver;
        solver.run(chain_, initcc_, times_);
        auto results = solver.results_;
        results.to_csv(resultpath_);
    }
    else
    {
        throw std::runtime_error(fmt::format("Invalid solver '{}'", solvertype_));
    }
}

std::string Model::chainpath() const
{
    return chainpath_;
}

void Model::set_chainpath(const std::string &chainpath)
{
    chainpath_ = chainpath;
    chain_ = Chain(chainpath_.c_str());
}
