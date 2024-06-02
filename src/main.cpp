#include "pugixml.hpp"
#include "model.h"
#include <fmt/core.h>

int main(int argc, char *argv[]) {
    (void) argc;

    std::string inputpath(argv[1]);
    fmt::print("Running {:s}...\n", inputpath);
    Model myinput(inputpath);
    myinput.run();

    return 0;
}
