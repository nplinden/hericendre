#include "input.h"
#include "chain.h"
#include "solver.h"
#include "utils.h"
#include "yaml-cpp/yaml.h"
#include <iostream>
#include <numeric>

Input::Input(std::string inputpath) {
  YAML::Node input =
      YAML::LoadFile("/home/nlinden/workspace/hericendre/data/example1.yaml");
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
  auto results = solver.run(chain, concentrations_, times_);

  for (const auto &N : results) {
    for (int i = 0; i < N.size(); i++) {
      if (N[i] > 1e-15)
        fmt::print("{:8s},{:.4e},{: .4e}\n", chain.nuclides_[i]->name_,
                   chain.nuclides_[i]->dconst_, N[i]);
    }
    double sum = std::accumulate(N.begin(), N.end(), 0.);
    fmt::print("sum = {:.4e}\n", sum);
  }
}
