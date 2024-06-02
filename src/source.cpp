#include <fmt/core.h>
#include <source.h>
#include <vector>

Source::Source(const pugi::xml_node &sourceNode) {
    type_ = sourceNode.attribute("type").value();
    particle_ = sourceNode.attribute("particle").value();
    interpolation_ = sourceNode.attribute("interpolation").value();

    auto params = std::string(sourceNode.child_value("parameters"));
    if (!params.empty()) {
        std::vector<std::string> params_vector;
        params_vector.emplace_back();
        for (const auto c: params) {
            if (c == ' ') {
                params_vector.emplace_back();
            } else {
                params_vector.back() += c;
            }
        }
        const size_t length = params_vector.size() / 2;
        std::vector<std::string> intensity(params_vector.begin() + length,
                                           params_vector.end());
        std::vector<std::string> energy(params_vector.begin(),
                                        length + params_vector.begin());
        for (auto &e: energy) {
            this->energy_.push_back(stod(e));
        }
        for (auto &i: intensity) {
            this->intensity_.push_back(stod(i));
        }
    }

    if (this->type_ == "mixture") {
        for (pugi::xml_node pairNode = sourceNode.child("pair"); pairNode;
             pairNode = pairNode.next_sibling("pair")) {
            std::string proba = pairNode.attribute("probability").value();
            pair_probabilities_.push_back(stod(proba));
            pairs_.emplace_back(pairNode.child("dist"));
        }
    }
}
