#ifndef NUCLIDE_HPP_INCLUDED
#define NUCLIDE_HPP_INCLUDED
#include "source.h"
#include <fmt/core.h>
#include <map>
#include <memory>
#include <pugixml.hpp>
#include <string>
#include <vector>
// #include "nfy.h"

class Decay;
class NReaction;
class Fission;
using SourcePtr = std::shared_ptr<Source>;

class Nuclide {
public:
 // CONSTRUCTORS
 Nuclide(const std::string &name, double dconst);

 /**
  * \brief Creates a Nuclide object using an nuclide-tagged xml node.
  *
  */
 explicit Nuclide(const pugi::xml_node &nuclideNode);

 // Members functions
 /**
  * \brief Add a nuclide node to an xml chain file.
  *
  * \param rootnode: The parent node to add to.
  */
 void addNode(pugi::xml_node &rootnode);

 /**
  * \brief Define an order on nuclide with a less than operator.
  */
 bool operator<(const Nuclide &other) const;

 /**
  * \brief Return the Nuclide object in string form.
  */
 std::string str() const;

 /**
  * \brief A static member function that return the the (z, a, m) tuple
  * from a nuclide name.
  *
  * \param std::string name: a nuclide name.
  */
 static std::tuple<int, int, int> getZam(const std::string &name);

 // MEMBER VARIABLES
 /**
  * \brief The vector of decays reaction that the nuclide can undergo.
  */
 std::vector<std::shared_ptr<Decay> > decays_;

 /**
  * \brief The vector of decays reaction that result in the nuclide's creation.
  */
 std::vector<std::shared_ptr<Decay> > decaysUp_;

 /**
  * \brief The vector of reactions that the nuclide can undergo.
  */
 std::vector<std::shared_ptr<NReaction> > reactions_;

 /**
  * \brief The vector of reactions that result in the nuclide's creation.
  */
 std::vector<std::shared_ptr<NReaction> > reactionsUp_;

 /**
  * \brief The vector of particles sources of the nuclide.
  */
 std::vector<SourcePtr> sources_;

 /**
  * \brief The decay constant of the nuclide.
  */
 double dconst_;
 double getDconst() const { return dconst_; };
 void setDconst(const double dconst) { dconst_ = dconst; };

 /**
  * \brief The decay energy of the nuclide.
  */
 double denergy_;
 double getDenergy() const { return denergy_; };

 /**
  * \brief The name of the nuclide.
  */
 std::string name_;
 std::string getName() const { return name_; };
 void setName(const std::string &name) { name_ = name; };

 // std::string nfyParent_ = "";

 /**
  * \brief The index of the isotope in a Chain object.
  */
 size_t idInChain;

 /**
  * \brief The number of reaction types that the nuclide can undergo.
  */
 int nreac_;

 /**
  * \brief The number of decay types that the nuclide can undergo.
  */
 int ndecay_;

 /**
  * \brief A mapping of element symbols to charge number.
  */
 static const std::map<std::string, int> ELEMENTS;

 /**
  * \brief The nuclide's charge number.
  */
 int z_;

 /**
  * \brief The nuclide's mass number.
  */
 int a_;

 /**
  * \brief The nuclide's metastable state.
  */
 int m_;

 /**
  * \brief The nuclide's id following the id = 10000 * Z + 10 * A + M
  * convention.
  */
 int zam_;
};

template<>
class fmt::formatter<Nuclide> {
public:
 constexpr auto parse(format_parse_context &ctx) { return ctx.begin(); }

 template<typename Context>
 constexpr auto format(Nuclide const &nuc, Context &ctx) const {
  return format_to(ctx.out(), "{}", nuc.name_);
 }
};
#endif
