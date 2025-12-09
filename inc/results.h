#ifndef RESULTS_HPP_INCLUDED
#define RESULTS_HPP_INCLUDED
#include <Eigen/Sparse>
#include <highfive/H5Easy.hpp>
#include <vector>

/**
 * @brief Container class for storing and exporting nuclide concentration results
 * 
 * This class holds the results of decay or depletion calculations, including
 * nuclide concentrations at different time points. It provides methods to export
 * the data to CSV and HDF5 formats.
 */
class Results {
public:
  std::vector<std::vector<double>> cc_; ///< Concentration matrix [time][nuclide]
  std::vector<std::string> nuclides_;  ///< Names of the nuclides
  std::vector<double> times_;          ///< Time points for the results

  /**
   * @brief Default constructor
   */
  Results();

  /**
   * @brief Construct Results from concentration vectors
   * 
   * @param cc Concentration matrix as nested vectors [time][nuclide]
   * @param nuclides Vector of nuclide names
   * @param times Vector of time points
   */
  Results(const std::vector<std::vector<double>> &cc,
          const std::vector<std::string> &nuclides,
          const std::vector<double> &times);

  /**
   * @brief Construct Results from Eigen vectors
   * 
   * @param cc Vector of Eigen concentration vectors, one per time point
   * @param nuclides Vector of nuclide names
   * @param times Vector of time points
   */
  Results(const std::vector<Eigen::VectorXd> &cc,
          const std::vector<std::string> &nuclides,
          const std::vector<double> &times);

  /**
   * @brief Export results to CSV file
   * 
   * Writes the concentration data to a CSV file with nuclides as columns
   * and time points as rows.
   * 
   * @param path Output file path
   * @param ignore_zeros If true, omit nuclides that have zero concentration at all times
   */
  void to_csv(const std::string &path, bool ignore_zeros = true) const;

  /**
   * @brief Export results to HDF5 file
   * 
   * Writes the concentration data, nuclide names, and time points to an HDF5 file.
   * 
   * @param file Open HDF5 file object to write to
   */
  void to_hdf5(H5Easy::File &file) const;
};

#endif
