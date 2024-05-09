#include "pugixml.hpp"
// #include "yaml-cpp/node/node.h"
#include "input.h"
#include <Eigen/Sparse>
#include <chain.h>
#include <decay.h>
#include <fmt/core.h>
#include <fmt/os.h>
#include <fmt/ranges.h>
#include <nuclide.h>

#include <solver.h>

int main(int argc, char *argv[]) {
  (void)argc;

  std::string inputpath(argv[1]);
  fmt::print("Running {:s}...\n", inputpath);
  Input myinput(inputpath);
  myinput.run();

  return 0;
}
