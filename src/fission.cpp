#include "fission.hpp"
#include "utils.hpp"
#include <fmt/core.h>
#include <fmt/os.h>
#include <fmt/ranges.h>

Fission::Fission(const pugi::xml_node &nfyNode, NuclidePtr parent) {
  parent_ = parent;
  parentName_ = parent->name_;

  this->energies_ = split(std::string(nfyNode.child_value("energies")));
  for (pugi::xml_node fyNode = nfyNode.child("fission_yields"); fyNode;
       fyNode = fyNode.next_sibling("fission_yields")) {
    std::string energy = fyNode.attribute("energy").value();
    fy_[energy] = std::map<std::string, double>();
    std::vector<std::string> products = split(fyNode.child_value("products"));
    std::vector<std::string> data = split(fyNode.child_value("data"));
    for (int i = 0; i < products.size(); i++) {
      fy_[energy][products[i]] = stod(data[i]);
    }
  }
}

void Fission::addNode(pugi::xml_node &nfyNode) {
  std::string energies_str = concatenate(this->energies_);
  auto energiesNode = nfyNode.append_child("energies");
  energiesNode.text() = energies_str.c_str();
  for (auto const &[e, v] : this->fy_) {
    auto fyNode = nfyNode.append_child("fission_yields");
    fyNode.append_attribute("energy") = e.c_str();
    std::string productStr = "";
    std::string dataStr = "";
    for (auto const &[target, targetfy] : v) {
      productStr += fmt::format("{} ", target);
      dataStr += fmt::format("{} ", targetfy);
    }
    productStr.pop_back();
    dataStr.pop_back();
    auto products = fyNode.append_child("products");
    products.text() = productStr.c_str();
    auto data = fyNode.append_child("data");
    data.text() = dataStr.c_str();
  }
}
