#include "model.h"
#include <fmt/core.h>
#include <highfive/highfive.hpp>

int main(int argc, char *argv[])
{
    fmt::print("    __  __          _                     __\n");
    fmt::print("   / / / /__  _____(_)_______  ____  ____/ /_______ \n");
    fmt::print("  / /_/ / _ \\/ ___/ / ___/ _ \\/ __ \\/ __  / ___/ _ \\\n");
    fmt::print(" / __  /  __/ /  / / /__/  __/ / / / /_/ / /  /  __/\n");
    fmt::print("/_/ /_/\\___/_/  /_/\\___/\\___/_/ /_/\\__,_/_/   \\___/\n");
    fmt::print("\n");

    if (argc < 2)
    {
        fmt::print("[ERROR] No input file was provided! Exiting.\n");
        return 0;
    }

    std::string inputPath(argv[1]);
    fmt::print("Running {:s}...\n", inputPath);
    Model inputModel(inputPath);
    inputModel.run();

    return 0;
}
