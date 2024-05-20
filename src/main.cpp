#include "pugixml.hpp"
// #include "yaml-cpp/node/node.h"
#include "decaysolver.h"
#include "input.h"
#include <Eigen/Sparse>
#include <chain.h>
#include <decay.h>
#include <fmt/core.h>
#include <fmt/os.h>
#include <fmt/ranges.h>
#include <map>
#include <nuclide.h>

#include <solver.h>

// int main(int argc, char *argv[]) {
//   (void)argc;
//
//   std::string inputpath(argv[1]);
//   fmt::print("Running {:s}...\n", inputpath);
//   Input myinput(inputpath);
//   myinput.run();
//
//   return 0;
// }

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;
  Chain chain("/home/nlinden/workspace/hericendre/data/chain-jeff33-long.xml");
  DecaySolver solver;

  std::map<std::string, double> initcc;
  for (const auto &nuclide : chain.nuclides_) {
    initcc[nuclide->name_] = 1.;
  }
  std::vector<double> times;
  for (double i = -10; i <= 17; i = i + 1)
    times.push_back(pow(10, i));

  solver.run(chain, initcc, times);
  solver.results_.to_csv("analytical_results.csv");
}
