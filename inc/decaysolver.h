#ifndef DECAYSOLVER_HPP_INCLUDED
#define DECAYSOLVER_HPP_INCLUDED
#include "chain.h"
#include "results.h"
#include <map>
#include <string>
#include <highfive/H5Easy.hpp>
#include <highfive/highfive.hpp>

/**
 * @class DecaySolver
 * @brief A class for solving decay chain problems.
 *
 * The `DecaySolver` class is responsible for simulating the behavior of a decay chain.
 * It computes coefficients for stable and unstable nuclides and solves the system
 * over a given set of time steps.
 */
class DecaySolver
{
public:
    /**
     * @brief Constructs a DecaySolver object.
     * @param chain A reference to the `Chain` object representing the decay chain.
     */
    explicit DecaySolver(Chain &chain);

    /**
     * @brief Runs the decay simulation over a series of time steps.
     * @param ccMap A map of nuclide names to their initial concentrations.
     * @param times A vector of time points at which to compute the solution.
     * @return A `Results` object containing the computed concentrations over time.
     */
    Results run(const std::map<std::string, double> &ccMap, std::vector<double> times);

private:
    /**
     * @brief Computes the coefficients for all nuclides in the decay chain.
     * @param ccMap A map of nuclide names to their initial concentrations.
     */
    void compute_coeffs(std::map<std::string, double> ccMap);

    /**
     * @brief Computes the `Fik` coefficients for a given nuclide.
     * @param nuclide A pointer to the nuclide for which to compute the coefficients.
     */
    void compute_Fik(const NuclidePtr nuclide);

    /**
     * @brief Computes the `Ns` coefficients for a given nuclide.
     * @param nuclide A pointer to the nuclide for which to compute the coefficients.
     * @param ccMap A map of nuclide names to their initial concentrations.
     */
    void compute_Ns(const NuclidePtr nuclide, const std::map<std::string, double> &ccMap);

    /**
     * @brief Computes the `Fii` coefficients for a given nuclide.
     * @param nuclide A pointer to the nuclide for which to compute the coefficients.
     * @param ccMap A map of nuclide names to their initial concentrations.
     */
    void compute_Fii(const NuclidePtr nuclide, const std::map<std::string, double> &ccMap);

    Chain chain_; ///< The decay chain being simulated.

    std::map<size_t, double> Ns;                  ///< Map of nuclide indices to their `Ns^{i}` coefficients.
    std::map<size_t, std::map<size_t, double>> F; ///< Nested map representing the `F_{i,j}` coefficients.
};

#endif
