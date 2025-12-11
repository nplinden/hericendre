#ifndef CHAIN_HPP_INCLUDED
#define CHAIN_HPP_INCLUDED
#include <Eigen/Sparse>
#include <decay.h>
#include <nuclide.h>
#include <string>
#include <vector>

using DecayPtr = std::shared_ptr<Decay>;
using NuclidePtr = std::shared_ptr<Nuclide>;

/**
 * @brief Represents a nuclear decay and transmutation chain.
 *
 * The Chain class manages a collection of nuclides and their interactions through
 * decay reactions, neutron-induced reactions, and fission processes. It provides
 * methods for querying nuclides, analyzing reachability, topological sorting, and
 * constructing decay matrices for depletion calculations.
 *
 * The chain can be loaded from OpenMC-format XML files or constructed programmatically.
 */
class Chain
{
public:
    // CONSTRUCTORS
    /**
     * @brief Constructs a chain from an OpenMC XML chain file.
     *
     * Parses the specified XML file to populate the chain with nuclides, decay modes,
     * neutron reactions, and fission yields.
     *
     * @param path Path to the depletion chain file in OpenMC XML format.
     */
    explicit Chain(const char *path);

    /**
     * @brief Constructs a chain from an OpenMC XML chain file.
     *
     * @param path Path to the depletion chain file in OpenMC XML format.
     */
    explicit Chain(const std::string &path);

    /**
     * @brief Constructs an empty chain.
     */
    Chain();

    // MEMBER FUNCTIONS
    /**
     * @brief Writes the chain to a file.
     *
     * @param path Path to the output file.
     * @return True if successful, false otherwise.
     */
    bool write(const char *path);

    // void restrict(std::vector<std::string> nuclides);
    // void removeNuclide(std::string nuc);

    /**
     * @brief Checks if a nuclide exists in the chain.
     *
     * @param name Name of the nuclide (e.g., "U235", "Xe135_m1").
     * @return True if the nuclide is present in the chain.
     */
    bool contains(const std::string &name) const;

    /**
     * @brief Finds a nuclide by its numeric identifier.
     *
     * @param nucid Nuclide ID calculated as: 10000*Z + 10*A + E,
     *              where Z is atomic number, A is mass number, E is isomeric state.
     * @return Shared pointer to the Nuclide object, or nullptr if not found.
     */
    NuclidePtr find(int nucid) const;

    /**
     * @brief Finds a nuclide by its name.
     *
     * @param name Name of the nuclide (e.g., "U235", "Xe135_m1").
     * @return Shared pointer to the Nuclide object, or nullptr if not found.
     */
    NuclidePtr find(const std::string &name) const;

    /**
     * @brief Finds the index of a nuclide by its numeric identifier.
     *
     * @param nucid Nuclide ID calculated as: 10000*Z + 10*A + E.
     * @return Index in the nuclides_ vector.
     * @throws std::out_of_range if nuclide not found.
     */
    size_t nuclide_index(int nucid) const;

    /**
     * @brief Finds the index of a nuclide by its name.
     *
     * @param name Name of the nuclide.
     * @return Index in the nuclides_ vector.
     * @throws std::out_of_range if nuclide not found.
     */
    size_t nuclide_index(const std::string &name) const;

    /**
     * @brief Saves the chain to a file.
     *
     * @param path Path to the output file.
     */
    void save(const std::string &path) const;

    /**
     * @brief Performs depth-first search traversal of the chain graph.
     *
     * Starting from a given nuclide, traverses all reachable nuclides through
     * decay and transmutation pathways, marking visited nodes.
     *
     * @param nucid Index of the starting nuclide in the nuclides_ vector.
     * @param visited Boolean vector tracking which nuclides have been visited.
     */
    void dfs(const size_t &nucid, std::vector<bool> &visited);

    /**
     * @brief Finds all nuclides reachable from a starting nuclide.
     *
     * Uses depth-first search to identify all nuclides that can be reached
     * through decay chains and neutron-induced transmutations.
     *
     * @param nucname Name of the starting nuclide.
     * @return Vector of names of all reachable nuclides.
     */
    std::vector<std::string> reachable(const std::string &nucname);

    /**
     * @brief Extracts a vector of all nuclide names in the chain.
     *
     * @return Vector of nuclide names in order of the nuclides_ vector.
     */
    std::vector<std::string> name_vector() const;

    /**
     * @brief Extracts a vector of all decay constants.
     *
     * @return Vector of decay constants (1/s) in order of the nuclides_ vector.
     */
    std::vector<double> dconst_vector() const;

    /**
     * @brief Sorts the chain in topological order.
     *
     * Reorders nuclides so that parent nuclides appear before their decay products,
     * which is beneficial for matrix solver performance.
     *
     * @return True if topological sorting was successful, false if cycles detected.
     */
    bool topological_sort();

    /**
     * @brief Perturbs duplicate decay constants to ensure uniqueness.
     *
     * Applies small perturbations (factor of 1E-14) to decay constants that have
     * identical values, which can cause numerical issues in some solver algorithms.
     */
    void tweak_dconst();

    /**
     * @brief Constructs the decay matrix for the Bateman equations.
     *
     * Builds a sparse matrix A where A[i,j] represents the rate at which nuclide j
     * produces nuclide i through decay. Diagonal elements are the negative of the
     * total decay constant. This matrix is used in dN/dt = A*N.
     *
     * @return Sparse matrix of dimension (n_nuclides × n_nuclides).
     */
    Eigen::SparseMatrix<double> decayMatrix() const;

    // MEMBER VARIABLES
    /**
     * @brief Collection of all nuclides in the chain.
     *
     * Each nuclide contains information about its properties, decay modes,
     * and neutron reaction cross sections. The index of each nuclide in this
     * vector is stored as idInChain for efficient lookup.
     */
    std::vector<NuclidePtr> nuclides_;

    /**
     * @brief Collection of all decay transitions in the chain.
     *
     * Each Decay object represents a radioactive decay mode (alpha, beta, etc.)
     * connecting a parent nuclide to a daughter nuclide with a specific branching ratio.
     */
    std::vector<DecayPtr> decays_;
};
#endif
