#include "model.h"
#include "chain.h"
#include "results.h"
#include "cramsolver.h"
#include "utils.h"
#include <cmath>
#include "yaml-cpp/yaml.h"
#include <fmt/ranges.h>
#include <algorithm>
#include <decaysolver.h>
#include <highfive/H5Easy.hpp>

Model::Model(const std::string &inputpath) {
    YAML::Node input = YAML::LoadFile(inputpath);
    chainpath_ = input["Chain"].as<std::string>();
    chain_ = Chain(chainpath_.c_str());
    resultpath_ = input["Results"].as<std::string>();

    readSolverType(input);
    readTimes(input);
    readCc(input);
}

Model::Model() = default;


void Model::readSolverType(const YAML::Node &input) {
    if (const auto modenode = input["Solver"]; !modenode) {
        solvertype_ = "CRAM48";
    } else {
        auto solvertype = modenode.as<std::string>();
        if (solvertype == "CRAM48") {
            solvertype_ = "CRAM48";
        } else if (solvertype == "Decay") {
            solvertype_ = "Decay";
        } else {
            throw std::runtime_error(fmt::format("Invalid solver type '{}'", solvertype));
        }
    }
}

void Model::readCc(const YAML::Node &input) {
    std::string ccmode = "Explicit";
    if (const auto mode = input["ConcentrationMode"]) ccmode = mode.as<std::string>();
    if (ccmode == "Explicit") {
        const std::vector<std::string> cc = split(input["Concentrations"].as<std::string>());
        for (size_t i = 0; i < cc.size(); i += 2)
            initcc_[cc[i]] = stod(cc[i + 1]);
    } else if (ccmode == "Uniform") {
        const double val = stod(input["Concentrations"].as<std::string>());
        for (const auto &nuclide: chain_.nuclides_)
            initcc_[nuclide->name_] = val;
    } else {
        throw std::runtime_error(fmt::format("Invalid ConcentrationMode '{}'", ccmode));
    }
}

void Model::readTimes(const YAML::Node &input) {
    const auto timemode = input["TimeMode"].as<std::string>();

    const std::vector<std::string> raw_times = split(input["Time"].as<std::string>());
    std::vector<std::vector<std::string> > lines;
    lines.emplace_back();
    double T = 0.;
    for (const auto &word: raw_times) {
        if (word == ";") {
            lines.emplace_back();
        } else {
            lines.back().push_back(word);
        }
    }

    for (const auto &line: lines) {
        if (line[0] == "linspace") {
            std::vector<double> values = linspace(line);
            times_.insert(times_.end(), values.begin(), values.end());
        } else if (line[0] == "logspace") {
            std::vector<double> values = logspace(line);
            times_.insert(times_.end(), values.begin(), values.end());
        } else {
            const double val = stod(line[0]);
            const double factor = (line.size() == 2) ? time_units.at(line[1]) : 1.;
            if (timemode == "Timestamps") {
                times_.push_back(val * factor);
            } else if (timemode == "Delta") {
                times_.push_back(T + val * factor);
                T = times_.back();
            }
        }
    }
    if (times_[0] != 0) {
        times_.insert(times_.begin(), 0.);
    }
    if (!std::is_sorted(times_.begin(), times_.end()))
        throw std::runtime_error("Time vector is not sorted");
}

std::vector<double> Model::linspace(const std::vector<std::string> &splat) const {
    const double start = stod(splat[1]);
    const double stop = stod(splat[2]);
    const int nstep = stoi(splat[3]);
    const std::string unit = (splat.size() == 5) ? splat[4] : "s";

    const double step = (stop - start) / nstep;
    std::vector<double> values;
    for (int istep = 0; istep < nstep; istep++) {
        values.push_back((start + istep * step) * time_units.at(unit));
    }
    return values;
}

std::vector<double> Model::logspace(const std::vector<std::string> &splat) const {
    const double start = std::log10(stod(splat[1]));
    const double stop = std::log10(stod(splat[2]));
    const int nstep = stoi(splat[3]);
    const std::string unit = (splat.size() == 5) ? splat[4] : "s";


    const double step = (stop - start) / nstep;
    std::vector<double> values;
    for (int istep = 0; istep < nstep; istep++) {
        values.push_back((pow(10, start + istep * step)) * time_units.at(unit));
    }
    return values;
}

void Model::run() {
    if (solvertype_ == "Decay") {
        DecaySolver solver(chain_);
        solver.run(initcc_, times_);
        auto results = solver.results_;
        results.to_csv(resultpath_);
        H5Easy::File file("results.h5", H5Easy::File::Overwrite);
        results.to_hdf5(file);
        solver.to_hdf5((file));
    } else if (solvertype_ == "CRAM48") {
        CRAMSolver solver;
        solver.run(chain_, initcc_, times_);
        auto results = solver.results_;
        results.to_csv(resultpath_);
    } else {
        throw std::runtime_error(fmt::format("Invalid solver '{}'", solvertype_));
    }
}

std::string Model::chainpath() const {
    return chainpath_;
}

void Model::set_chainpath(const std::string &chainpath) {
    chainpath_ = chainpath;
    chain_ = Chain(chainpath_.c_str());
}
