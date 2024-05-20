#include "chain.h"
#include "pugixml.hpp"
#include <Eigen/Sparse>
#include <fmt/core.h>
#include <fmt/os.h>
#include <iostream>

Chain::Chain(const char *path) {
  fmt::print("entering Chain\n");
  pugi::xml_document doc;
  doc.load_file(path);
  pugi::xml_node chainxml = doc.child("depletion_chain");

  // INTIALIZATION STEP
  std::cout << "Initialization:\n";
  for (pugi::xml_node nuclide = chainxml.child("nuclide"); nuclide;
       nuclide = nuclide.next_sibling("nuclide")) {

    nuclides_.push_back(std::make_shared<Nuclide>(Nuclide(nuclide)));

    nuclides_.back()->idInChain = nuclides_.size() - 1;
    for (pugi::xml_node decayNode = nuclide.child("decay"); decayNode;
         decayNode = decayNode.next_sibling("decay")) {

      decays_.push_back(
          std::make_shared<Decay>(Decay(decayNode, nuclides_.back())));

      std::string type = decays_.back()->type_;
      double br = decays_.back()->branchingRatio_;
      if (Decay::SECONDARIES.find(type) != Decay::SECONDARIES.end()) {
        auto secVector = Decay::SECONDARIES.at(type);
        std::map<std::string, double> multiplicities;
        for (const auto &target : secVector) {
          if (multiplicities.find(target) != multiplicities.end()) {
            multiplicities[target] += br;
          } else {
            multiplicities[target] = br;
          }
        }

        for (const auto &kv : multiplicities) {
          fmt::print("{} -> [{} {}] -> {}\n", nuclides_.back()->name_, type,
                     kv.second, kv.first);
          decays_.push_back(std::make_shared<Decay>(
              Decay(type, kv.first, kv.second, nuclides_.back())));
        }
      }
    }
    // for (pugi::xml_node reactionNode = nuclide.child("reaction");
    // reactionNode;
    //      reactionNode = reactionNode.next_sibling("reaction")) {
    //   if (reactionNode.attribute("type").value() != std::string("fission")) {
    //     reactions_.push_back(std::make_shared<NReaction>(
    //         NReaction(reactionNode, nuclides_.back())));
    //   }
    // }

    // for (pugi::xml_node nfyNode = nuclide.child("neutron_fission_yields");
    //      nfyNode; nfyNode = nfyNode.next_sibling("neutron_fission_yields")) {
    //   fissions_.push_back(
    //       std::make_shared<Fission>(Fission(nfyNode, nuclides_.back())));
    // }
  }
  std::cout << "\t" << nuclides_.size() << " nuclides found\n";
  std::cout << "\t" << decays_.size() << " decay reactions found\n";
  std::cout << "\t" << reactions_.size() << " neutron reactions found\n";
  // std::cout << "\t" << fissions_.size() << " fissions found\n";

  std::cout << "Binding decays\n";
  // BINDING STEP
  for (auto &dec : decays_) {
    dec->parent_->decays_.push_back(dec);
    if (!dec->targetName_.empty()) {
      dec->target_ = this->find(dec->targetName_);
      dec->target_->decaysUp_.push_back(dec);
    }
  }

  // std::cout << "Binding reaction\n";
  // for (auto &reac : reactions_) {
  //   reac->parent_->reactions_.push_back(reac);
  //   if (!reac->targetName_.empty())
  //     reac->target_ = this->find(reac->targetName_);
  //   reac->target_->reactionsUp_.push_back(reac);
  // }

  // for (auto& nuc : nuclides_){
  //   if (!nuc->nfyParent_.empty()){
  //     NuclidePtr target = this->find(nuc->nfyParent_) ;
  //     nuc->neutronFissionYields_ = target->neutronFissionYields_ ;
  //   }
  // }

  // std::cout << "Check...";
  // // COHERENCE CHECK
  // for (auto &n : nuclides_) {
  //   if (n->ndecay_ != n->decays_.size())
  //     throw std::invalid_argument(
  //         fmt::format("{} {} {}", n->name_, n->ndecay_, n->decays_.size()));
  //   if (n->nreac_ != n->reactions_.size())
  //     throw std::invalid_argument(
  //         fmt::format("{} {} {}", n->name_, n->nreac_,
  //         n->reactions_.size()));
  // }
  // std::cout << "done" << std::endl;
}

Chain::Chain() {}

void Chain::write(const char *path) {
  pugi::xml_document doc;
  auto root = doc.append_child("depletion_chain");

  for (auto nuclide : this->nuclides_) {
    nuclide->addNode(root);
  }
  doc.save_file(path, "  ");
}

NuclidePtr Chain::find(int nucid) const {
  for (auto nuc : nuclides_) {
    if (nuc->zam_ == nucid) {
      return nuc;
    }
  }
  std::string err_msg =
      fmt::format("Nuclide zam {} does not exist in the chain", nucid);
  throw std::invalid_argument(err_msg);
}

NuclidePtr Chain::find(std::string name) const {
  for (auto nuc : nuclides_) {
    if (nuc->name_ == name) {
      return nuc;
    }
  }
  std::string err_msg = fmt::format(
      "[find(std::string name)] Nuclide {} does not exist in the chain", name);
  throw std::invalid_argument(err_msg);
}

int Chain::nuclide_index(int nucid) const {
  for (size_t i = 0; i < nuclides_.size(); i++) {
    if (nuclides_[i]->zam_ == nucid)
      return i;
  }
  std::string err_msg =
      fmt::format("Nuclide zam {} does not exist in the chain", nucid);
  throw std::invalid_argument(err_msg);
}

int Chain::nuclide_index(std::string name) const {
  for (size_t i = 0; i < nuclides_.size(); i++) {
    if (nuclides_[i]->name_ == name)
      return i;
  }
  std::string err_msg = fmt::format("[nuclide_index(std::string name)] Nuclide "
                                    "{} does not exist in the chain",
                                    name);
  throw std::invalid_argument(err_msg);
}

Eigen::SparseMatrix<double> Chain::decayMatrix() const {
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(10000);
  for (size_t inuc = 0; inuc < this->nuclides_.size(); inuc++) {
    NuclidePtr nuc = this->nuclides_[inuc];
    triplets.push_back(Eigen::Triplet<double>(inuc, inuc, -nuc->dconst_));
    for (const auto &d : nuc->decays_) {
      if (!d->targetName_.empty()) {
        int jnuc = this->nuclide_index(d->targetName_);
        triplets.push_back(Eigen::Triplet<double>(
            jnuc, inuc, nuc->dconst_ * d->branchingRatio_));
      }
    }
  }
  Eigen::SparseMatrix<double> M(this->nuclides_.size(), this->nuclides_.size());
  M.setFromTriplets(triplets.begin(), triplets.end());
  M.makeCompressed();
  return M;
}

void Chain::dfs(int nucid, std::vector<bool> &visited) {
  if (visited[nucid])
    return;
  visited[nucid] = true;

  std::vector<DecayPtr> decays = nuclides_[nucid]->decays_;
  std::vector<int> neighbours;
  // fmt::print("{}: {}\n", nucid, nuclides_[nucid]->name_);
  for (DecayPtr decay : nuclides_[nucid]->decays_) {
    // fmt::print("\t{}\n", decay->target_->name_);
    neighbours.push_back(decay->target_->idInChain);
  }
  for (const auto neighboursId : neighbours) {
    this->dfs(neighboursId, visited);
  }
}

std::vector<std::string> Chain::reachable(std::string nucname) {
  int n = nuclides_.size();
  std::vector<bool> visited;
  for (int i = 0; i < n; i++)
    visited.push_back(false);
  int initialId = this->find(nucname)->idInChain;
  this->dfs(initialId, visited);

  std::vector<std::string> names;
  for (int i = 0; i < n; i++) {
    if (visited[i])
      names.push_back(nuclides_[i]->name_);
  }
  return names;
}

void Chain::dump_matrix(std::string path) {
  Eigen::MatrixXd dMat(this->decayMatrix());
  fmt::print("size = [{}, {}]\n", dMat.rows(), dMat.cols());

  auto out = fmt::output_file(path);
  out.print(",");
  for (size_t inuc = 0; inuc < this->nuclides_.size(); inuc++) {
    out.print("{}", this->nuclides_[inuc]->name_);
    if (inuc < this->nuclides_.size() - 1)
      out.print(",");
  }
  out.print("\n");
  for (long int row = 0; row < dMat.rows(); row++) {
    out.print("{},", this->nuclides_[row]->name_);
    for (long int col = 0; col < dMat.cols(); col++) {
      out.print("{}", dMat(row, col));
      if (col < dMat.cols() - 1)
        out.print(",");
    }
    out.print("\n");
  }
  out.close();
}

// void Chain::removeNuclide(std::string nuc) {}

// void Chain::restrict(std::vector<std::string> allowed) {
// std::vector<NuclidePtr> newNuclides ;
// for (auto& nuc : nuclides_){
//   if (std::find(allowed.begin(), allowed.end(), nuc->name_) !=
//   allowed.end()){
//     newNuclides.push_back(nuc) ;
//   }
// }

// std::vector<DecayPtr> newDecays ;
// for (auto& decay : decays_){
//   bool parentOk = std::find(allowed.begin(), allowed.end(),
//   decay->parentName_) != allowed.end() ; bool targetOk =
//   std::find(allowed.begin(), allowed.end(), decay->targetName_) !=
//   allowed.end() ; if (parentOk && targetOk){
//     newDecays.push_back(decay) ;
//   }
// }
// }
