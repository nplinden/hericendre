#include "model.h"
#include <fmt/core.h>
#include <highfive/highfive.hpp>
#include <microxs.h>

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
        return EXIT_FAILURE;
    }

    std::string inputPath(argv[1]);

    Model model;
    try
    {
        model = Model(inputPath);
        model.summarize();
    }
    catch (const std::exception &e)
    {
        fmt::print(stderr, "[ERROR] Failed to initialize model: {}\n", e.what());
        return EXIT_FAILURE;
    }
    fmt::print("{}", model.microxs_.getXS("U235", "fission"));

    auto M = model.chain_.DepletionMatrix(model.microxs_, 1.);
    return EXIT_SUCCESS;

    try
    {
        model.run();
    }
    catch (const std::exception &e)
    {
        fmt::print(stderr, "[ERROR] Simulation run failed: {}\n", e.what());
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
