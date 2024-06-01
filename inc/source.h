#ifndef SOURCE_HPP_INCLUDED
#define SOURCE_HPP_INCLUDED
#include <pugixml.hpp>
#include <string>
#include <vector>

class Source {
public:
  /**
   * \brief Creates a Source object using an source-tagged xml node.
   *
   */
  explicit Source(const pugi::xml_node &sourceNode);
  /**
   * \brief The type of source, "discrete", "tabulated", etc.
   */
  std::string type_;

  /**
   * \brief The emitted particle type.
   */
  std::string particle_;

  /**
   * \brief The interpolation scheme to be used.
   */
  std::string interpolation_;

  /**
   * \brief the energy values for which the source is defined.
   */
  std::vector<double> energy_;

  /**
   * \brief Intensity of the source, should be of same length as the energy_
   * attribute.
   */
  std::vector<double> intensity_;

  /**
   * \brief For pair source representation, this contains the relevent Source
   * object.
   */
  std::vector<Source> pairs_;

  /**
   * \brief For pair source representation, this contains the probabilities for
   * each source type.
   */
  std::vector<double> pair_probabilities_;
};

#endif
