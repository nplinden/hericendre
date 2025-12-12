#include "model.h"
#include "chain.h"
#include "results.h"
#include "cramsolver.h"
#include "utils.h"
#include <cmath>
#include <fmt/ranges.h>
#include <fmt/core.h>
#include <algorithm>
#include <decaysolver.h>
#include <highfive/H5Easy.hpp>
#include <toml++/toml.hpp>
#include <optional>
#include <set>
#include <filesystem>

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
    readTime(tbl);
    readMaterial(tbl);
}

Model::Model() = default;

void Model::readSettings(const toml::table &tbl)
{
    // DEPLETION CHAIN
    std::optional<std::string> chainpath = tbl["Settings"]["chain"].value<std::string>();
    if (chainpath)
    {
        chainpath_ = *chainpath;
        const std::filesystem::path p{chainpath_};
        if (!std::filesystem::exists(p))
        {
            std::string msg = fmt::format("Depletion chain file {} does not exist", chainpath_);
            throw std::runtime_error(msg);
        }
        this->chain_ = Chain(chainpath_.c_str());
    }
    else
    {
        throw std::runtime_error("No depletion chain file was specified");
    }

    // RESULT PATH
    std::optional<std::string> resultpath = tbl["Settings"]["results"].value<std::string>();
    if (resultpath)
    {
        this->resultpath_ = *resultpath;
    }
    else
    {
        throw std::runtime_error("No result file was specified");
    }

    std::optional<std::string> save_matrix = tbl["Settings"]["save_matrix"].value<std::string>();
    if (save_matrix)
    {
        this->chain_.save(*save_matrix);
    }

    // SOLVER TYPE
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
        std::string msg = fmt::format("Invalid solver type \"{}\", allowed solver: {}", *solver, fmt::join(allowed_solvers, ", "));
        throw std::runtime_error(msg);
    }
}

void Model::readTime(const toml::table &tbl)
{
    // TIMESTAMPS
    auto timestamps = tbl["Time"]["timestamps"].as_array();
    if (!timestamps)
    {
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
            throw std::runtime_error("Invalid time unit name.");
        }
        else
        {
            this->time_unit_name_ = *name.value<std::string>();
        }
        auto magnitude = unit_table["magnitude"];
        if (!magnitude)
        {
            throw std::runtime_error("Invalid time unit magnitude.");
        }
        else if (magnitude.is<double>() || magnitude.is<int64_t>())
        {
            this->time_unit_magnitude_ = *magnitude.value<double>();
        }
        else
        {
            throw std::runtime_error("Invalid time unit magnitude.");
        }
    }
    else
    {
        throw std::runtime_error("Time unit must be either a string or a table.");
    }

    for (auto &t : this->times_)
    {
        t *= this->time_unit_magnitude_;
    }

    if (!std::is_sorted(this->times_.begin(), this->times_.end()))
    {
        throw std::runtime_error("Time vector is not sorted.");
    }
}

void Model::readMaterial(const toml::table &tbl)
{
    std::optional<std::string> xspath = tbl["Material"]["microxs"].value<std::string>();
    if (xspath)
        microxs_ = MicroXS(*xspath);

    std::optional<double> uniform = tbl["Material"]["uniform"].value<double>();
    auto concentrations = tbl["Material"]["concentrations"].as_table();
    if (concentrations && uniform)
    {
        throw std::runtime_error("uniform and concentration can't be defined at the same time");
    }
    else if (uniform)
    {
        for (const auto &nuclide : chain_.nuclides_)
            initcc_[nuclide->name_] = *uniform;
    }
    else if (concentrations)
    {
        for (const auto &nuclide : *concentrations)
        {
            auto conc = nuclide.second.value<double>();
            if (!conc)
            {
                throw std::runtime_error("Concentration must be a number");
            }
            else
            {
                initcc_[std::string(nuclide.first.str())] = *conc;
            }
        }
    }
}

std::vector<double> Model::compute_time_function(const std::string &str)
{
    std::vector<std::string> splat = split(str);
    if (splat.size() < 4)
    {
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
            throw std::runtime_error("Invalid time function in timestamps definition");
        }
    }
}

std::vector<double> Model::linspace(const std::vector<std::string> &splat) const
{
    const double start = stod(splat[1]);
    const double stop = stod(splat[2]);
    const int nstep = stoi(splat[3]) - 1;

    const double step = (stop - start) / nstep;
    std::vector<double> values;
    for (int istep = 0; istep <= nstep; istep++)
    {
        values.push_back((start + istep * step));
    }
    return values;
}

std::vector<double> Model::logspace(const std::vector<std::string> &splat) const
{
    const double start = stod(splat[1]);
    const double stop = stod(splat[2]);
    const int nstep = stoi(splat[3]) - 1;

    const double step = (stop - start) / nstep;
    std::vector<double> values;
    for (int istep = 0; istep <= nstep; istep++)
    {
        values.push_back((pow(10, start + istep * step)));
    }
    return values;
}

void Model::run()
{
    Results results;
    if (solvertype_ == "Decay")
    {
        DecaySolver solver(chain_);
        results = solver.run(initcc_, times_);
    }
    else if (solvertype_ == "CRAM48")
    {
        CRAMSolver solver;
        results = solver.run(chain_, initcc_, times_);
    }
    else
    {
        throw std::runtime_error(fmt::format("Invalid solver '{}'", solvertype_));
    }
    auto extension = std::filesystem::path(this->resultpath_).extension();
    if (extension == ".h5" || extension == ".hdf5")
    {
        H5Easy::File file(this->resultpath_, H5Easy::File::Overwrite);
        results.to_hdf5(file);
    }
    else
    {
        results.to_csv(this->resultpath_);
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

void Model::summarize()
{
    int width = 50;
    fmt::print("┌{0:─^{1}}┐\n", this->name_, width); // Title
    fmt::print("│{1:.<{0}}{2:.>{0}}│\n", width / 2, "Name", this->name_);

    if (this->chainpath_.size() > 40)
    {
        auto filename = std::filesystem::path(this->chainpath_).filename().string();
        fmt::print("│{1:.<{0}}{2:.>{0}}│\n", width / 2, "Chain", filename);
    }
    else
    {
        fmt::print("│{1:.<{0}}{2:.>{0}}│\n", width / 2, "Chain", this->chainpath_);
    }
    fmt::print("│{1:.<{0}}{2:.>{0}}│\n", width / 2, "Number of nuclides", this->chain_.nuclides_.size());
    fmt::print("│{1:.<{0}}{2:.>{0}}│\n", width / 2, "Solver", this->solvertype_);

    fmt::print("│{1:.<{0}}{2:.>{0}}│\n", width / 2, "Number of timestamps", this->times_.size());
    fmt::print("│{1:.<{0}}{2:.>{0}}│\n", width / 2, "Begin date", this->times_.front());
    fmt::print("│{1:.<{0}}{2:.>{0}}│\n", width / 2, "End date", this->times_.back());

    fmt::print("└{0:─^{1}}┘\n", "", width);
}
