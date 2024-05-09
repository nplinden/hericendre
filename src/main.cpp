#include "pugixml.hpp"
// #include "yaml-cpp/node/node.h"
#include "input.h"
#include "yaml-cpp/yaml.h"
#include <Eigen/Sparse>
#include <chain.h>
#include <decay.h>
#include <fmt/core.h>
#include <fmt/os.h>
#include <fmt/ranges.h>
#include <highfive/H5Easy.hpp>
#include <iostream>
#include <map>
#include <nuclide.h>
#include <numeric>

#include <solver.h>
#include <vector>

int main(int argc, char *argv[]) {
  std::string inputpath(argv[1]) ;
  fmt::print("Running {:s}...\n", inputpath);
  Input myinput(inputpath);
  // fmt::print("powerlevel = {}\n", myinput.powerlevel_);
  // fmt::print("times = {}\n", myinput.times_);
  // fmt::print("cc = {}\n", myinput.concentrations_);
  myinput.run();

  return 0;
}
