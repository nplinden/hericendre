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
  * \param path: path to a depletion chain.
  */
 explicit Chain(const char *path);

 explicit Chain(const std::string &path);

 Chain();

 // MEMBER FUNCTIONS
 bool write(const char *path);

 // void restrict(std::vector<std::string> nuclides);
 // void removeNuclide(std::string nuc);

 bool isIn(const std::string &name) const;

 /**
  * \brief Returns a pointer to the desired Nuclide object, as defined by
  * its id.
  *
  * \param nucid: the id of the nuclide. id = 10000 * Z + 10 * A + E
  */
 NuclidePtr find(int nucid) const;

 /**
  * \brief Returns a pointer to the desired Nuclide object, as defined by
  * its name.
  *
  * \param name: the name of the nuclide.
  */
 NuclidePtr find(const std::string &name) const;

 /**
  * \brief Finds the index of a nuclide as defined by its id in the nuclide
  * vector.
  *
  * \param nucid: the id of the nuclide. id = 10000 * Z + 10 * A + E
  */
 size_t nuclide_index(int nucid) const;

 /**
  * \brief Finds the index of a nuclide as defined by its name in the nuclide
  * vector.
  *
  * \param name: the name of the nuclide.
  */
 size_t nuclide_index(const std::string &name) const;

 /**
  * \brief Dump the decay matrix to a csv file.
  *
  * \param path: Path of the file to write the matrix in.
  */
 void dump_matrix(const std::string &path) const;

 /**
  * \brief An implementation of depth first search algorithm to find all
  * visited nuclide nodes starting from an initial node.
  *
  * \param nucid: The id of the initial nuclide.
  * \param visited: A node visitation record.
  */
 void dfs(const size_t &nucid, std::vector<bool> &visited);

 /**
  * \brief Returns the list of all reachable nuclides names.
  *
  * \param nucname: the name of the nuclide.
  */
 std::vector<std::string> reachable(const std::string &nucname);

 std::vector<std::string> name_vector() const;
 std::vector<double> dconst_vector() const;

 /**
  * \brief Sort the chain in the topological order if possible.
  *
  */
 bool topological_sort();

 /**
  * \brief Tweaks decay constants of equal values by a factor of 1E-14.
  *
  */
 void tweak_dconst();

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
