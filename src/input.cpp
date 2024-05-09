#include "input.h"
#include "chain.h"
#include "results.h"
#include "solver.h"
#include "utils.h"
#include "yaml-cpp/yaml.h"
#include <iostream>
#include <numeric>

Input::Input(std::string inputpath) {
  YAML::Node input = YAML::LoadFile(inputpath);
  chainpath_ = input["chain"].as<std::string>();
  cyclemode_ = input["cycle_mode"].as<std::string>();
  std::vector<std::string> cycle = split(input["cycle"].as<std::string>());
  std::vector<double> times;
  for (int i = 0; i < cycle.size(); i++) {
    if (i % 2 == 0)
      times.push_back(std::stod(cycle[i]));
    else
      powerlevel_.push_back(std::stod(cycle[i]));
  }
  if (cyclemode_ == "t") {
    times_ = times;
  } else if (cyclemode_ == "dt") {
    times_.push_back(0.);
    for (auto &dt : times) {
      times_.push_back(times_.back() + dt);
    }
  };

  std::vector<std::string> cc =
      split(input["concentrations"].as<std::string>());

  for (int i = 0; i < cc.size(); i += 2) {
    concentrations_[cc[i]] = stod(cc[i + 1]);
  }
}

void Input::run() {
  Chain chain(chainpath_.c_str());
  Solver solver;
  solver.run(chain, concentrations_, times_);
  auto results = solver.results_;
  results.to_csv("results.csv");
}
