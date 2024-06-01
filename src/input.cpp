#include "input.h"
#include "chain.h"
#include "results.h"
#include "solver.h"
#include "utils.h"
#include <cmath>
#include "yaml-cpp/yaml.h"
#include <fmt/ranges.h>
#include <algorithm>

Input::Input(const std::string &inputpath) {
  YAML::Node input = YAML::LoadFile(inputpath);
  chainpath_ = input["Solver"].as<std::string>();
  chainpath_ = input["Chain"].as<std::string>();
  const auto timemode = input["TimeMode"].as<std::string>();
  readTimes(input, timemode);

  const std::vector<std::string> cc =
      split(input["Concentrations"].as<std::string>());

  for (size_t i = 0; i < cc.size(); i += 2) {
    concentrations_[cc[i]] = stod(cc[i + 1]);
  }
}

void Input::readTimes(const YAML::Node& input, const std::string& timemode) {
  const std::vector<std::string> raw_times = split(input["Time"].as<std::string>());
  std::vector<std::vector<std::string>> lines;
  lines.emplace_back();
  double T = 0. ;
  for (const auto& word: raw_times) {
    if (word == ";") {
      lines.emplace_back() ;
    } else {
      lines.back().push_back(word) ;
    }
  }

  for (const auto& line: lines) {
    if (line[0] == "linspace") {
      std::vector<double> values = linspace(line) ;
      times_.insert(times_.end(), values.begin(), values.end()) ;
    } else if (line[0] == "logspace") {
      std::vector<double> values = logspace(line) ;
      times_.insert(times_.end(), values.begin(), values.end()) ;
    } else {
      const double val = stod(line[0]);
      const double factor = (line.size() == 2) ? time_units.at(line[1]): 1.;
      if (timemode == "Timestamps") {
        times_.push_back(val * factor);
      } else if (timemode == "Delta") {
        times_.push_back(T + val * factor);
        T = times_.back() ;
      }
    }
  }
  if (times_[0] != 0) {
    times_.insert(times_.begin(), 0.) ;
  }
  if (!std::is_sorted(times_.begin(), times_.end()))
    throw std::runtime_error("Time vector is not sorted") ;

}

std::vector<double> Input::linspace(const std::vector<std::string>& splat) const {
  const double start = stod(splat[1]);
  const double stop = stod(splat[2]);
  const int nstep = stoi(splat[3]);
  const std::string unit = (splat.size() == 5) ? splat[4]: "s";

  const double step = (stop - start) / nstep;
  std::vector<double> values;
  for (int istep = 0; istep < nstep; istep++) {
    values.push_back((start + istep * step) * time_units.at(unit));
  }
  return values;
}

 std::vector<double> Input::logspace(const std::vector<std::string>& splat) const {
  const double start = std::log10(stod(splat[1]));
  const double stop = std::log10(stod(splat[2]));
  const int nstep = stoi(splat[3]);
  const std::string unit = (splat.size() == 5) ? splat[4]: "s";


  const double step = (stop - start) / nstep;
  std::vector<double> values;
  for (int istep = 0; istep < nstep; istep++) {
    values.push_back((pow(10, start + istep * step)) * time_units.at(unit));
  }
  return values;
}


void Input::run() const {
  const Chain chain(chainpath_.c_str());
  Solver solver;
  solver.run(chain, concentrations_, times_);
  auto results = solver.results_;
  results.to_csv("results.csv");
}
