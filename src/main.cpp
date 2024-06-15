#include "model.h"
#include <fmt/core.h>
#include <highfive/highfive.hpp>

int main(int argc, char *argv[]) {
    (void) argc;

    std::string inputpath(argv[1]);
    fmt::print("Running {:s}...\n", inputpath);
    Model myinput(inputpath);
    myinput.run();

    auto file = HighFive::File("foo.h5", HighFive::File::Create);
    return 0;
}
