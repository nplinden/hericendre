#ifndef CHAIN_HPP_INCLUDED
#define CHAIN_HPP_INCLUDED
#include <Eigen/Sparse>
#include <decay.h>
#include <fission.h>
#include <nreaction.h>
#include <nuclide.h>
#include <string>
#include <vector>

using DecayPtr = std::shared_ptr<Decay>;
using NReactionPtr = std::shared_ptr<NReaction>;
using NuclidePtr = std::shared_ptr<Nuclide>;
using NFYPtr = std::shared_ptr<Fission>;
class Chain {
public:
  // CONSTRUCTORS
  /**
   * \brief Builds a chain object from a path to a chain file in the OpenMC xml
   * format.
   *
   * \param int nucid: the id of the nuclide. id = 10000 * Z + 10 * A + E
   */
  Chain(const char *path);
  Chain();

  // MEMBER FUNCTIONS
  void write(const char *path);

  // void restrict(std::vector<std::string> nuclides);
  // void removeNuclide(std::string nuc);

  /**
   * \brief Returns a pointer to the desired Nuclide object, as defined by
   * its id.
   *
   * \param int nucid: the id of the nuclide. id = 10000 * Z + 10 * A + E
   */
  NuclidePtr find(int nucid) const;

  /**
   * \brief Returns a pointer to the desired Nuclide object, as defined by
   * its name.
   *
   * \param std::string name: the name of the nuclide.
   */
  NuclidePtr find(std::string name) const;

  /**
   * \brief Finds the index of a nuclide as defined by its id in the nuclide
   * vector.
   *
   * \param int nucid: the id of the nuclide. id = 10000 * Z + 10 * A + E
   */
  int nuclide_index(int nucid) const;

  /**
   * \brief Finds the index of a nuclide as defined by its name in the nuclide
   * vector.
   *
   * \param std::string name: the name of the nuclide.
   */
  int nuclide_index(std::string name) const;

  /**
   * \brief An implementation of depth first search algorithme to find all
   * visited nuclide nodes starting from an initial node.
   *
   * \param nucid: The id of the initial nuclide.
   * \param visited: A node visitation record.
   */
  void dfs(int nucid, std::vector<bool> &visited);

  /**
   * \brief Returns the list of all reachable nuclides names.
   *
   * \param std::string name: the name of the nuclide.
   */
  std::vector<std::string> reachable(std::string nucname);

  Eigen::SparseMatrix<double> decayMatrix() const;

  // MEMBER VARIABLES
  /**
   * \brief A vector of pointers to Nuclide objects. This constitutes the
   * list of nuclides available in the chain.
   */
  std::vector<NuclidePtr> nuclides_;

  /**
   * \brief A vector of pointers to Decay object. This consitutes the list
   * of decay reactions available in the chain.
   */
  std::vector<DecayPtr> decays_;

  /**
   * \brief A vector of pointers to NReaction object. This consitutes the list
   * of neutron reactions available in the chain.
   */
  std::vector<NReactionPtr> reactions_;

  /**
   * \brief A vector of pointers to Fission object. This consitutes the list
   * of fission reactions available in the chain.
   */
  std::vector<NFYPtr> fissions_;
};

#endif
