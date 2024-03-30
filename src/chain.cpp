#include <fmt/core.h>
#include <fmt/os.h>
#include <algorithm>
#include "chain.hpp"
#include "pugixml.hpp"
#include <iostream>
#include <cassert>
#include <Eigen/Sparse>
#include "source.hpp"

Chain::Chain(const char* path){
  pugi::xml_document doc ;
  pugi::xml_parse_result result = doc.load_file(path) ;
  pugi::xml_node chainxml = doc.child("depletion_chain") ;
  
  // INTIALIZATION STEP
  for (pugi::xml_node nuclide = chainxml.child("nuclide"); nuclide; nuclide = nuclide.next_sibling("nuclide"))
  {
    nuclides_.push_back(std::make_shared<Nuclide>(Nuclide(nuclide))) ;
    nuclides_.back()->idInChain = nuclides_.size() - 1 ;
    for (pugi::xml_node decayNode = nuclide.child("decay"); decayNode; decayNode = decayNode.next_sibling("decay")){
      decays_.push_back(std::make_shared<Decay>(Decay(decayNode, nuclides_.back()))) ;
    }
    for (pugi::xml_node reactionNode = nuclide.child("reaction"); reactionNode; reactionNode = reactionNode.next_sibling("reaction")){
      reactions_.push_back(std::make_shared<NReaction>(NReaction(reactionNode, nuclides_.back()))) ;
    }
    for (pugi::xml_node nfyNode = nuclide.child("neutron_fission_yields"); nfyNode; nfyNode = nfyNode.next_sibling("neutron_fission_yields")){
      fissions_.push_back(std::make_shared<Fission>(Fission(nfyNode, nuclides_.back()))) ;
    }
  }

  // BINDING STEP
  for (auto& dec : decays_){
    dec->parent_->decays_.push_back(dec) ;
    if (!dec->targetName_.empty())
      dec->target_ = this->find(dec->targetName_) ;
    dec->target_->decaysUp_.push_back(dec) ;
  }

  for (auto& reac : reactions_){
    reac->parent_->reactions_.push_back(reac) ;
    if (!reac->targetName_.empty())
      reac->target_ = this->find(reac->targetName_) ;
    reac->target_->reactionsUp_.push_back(reac) ;
  }

  

  // for (auto& nuc : nuclides_){
  //   if (!nuc->nfyParent_.empty()){
  //     NuclidePtr target = this->find(nuc->nfyParent_) ;
  //     nuc->neutronFissionYields_ = target->neutronFissionYields_ ;
  //   }
  // }

  // COHERENCE CHECK
  for (auto& n : nuclides_){
    if (n->ndecay_ != n->decays_.size())
      throw std::invalid_argument(fmt::format("{} {} {}", n->name_, n->ndecay_, n->decays_.size())) ;
    if (n->nreac_ != n->reactions_.size())
      throw std::invalid_argument(fmt::format("{} {} {}", n->name_, n->nreac_, n->reactions_.size())) ;
  }
}

Chain::Chain(){
  
}

void Chain::write(const char* path){
  pugi::xml_document doc;
  auto root = doc.append_child("depletion_chain");

  for (auto nuclide : this->nuclides_) {
    nuclide->addNode(root) ;
  }
  doc.save_file(path, "  ") ;
}

NuclidePtr Chain::find(int nucid) const{
    for (auto nuc : nuclides_){
      if (nuc->zam_ == nucid){
        return nuc ;
      }
    }
    std::string err_msg = fmt::format("Nuclide zam {} does not exist in the chain", nucid) ;
    throw std::invalid_argument(err_msg) ;
}

NuclidePtr Chain::find(std::string name) const{
    for (auto nuc : nuclides_){
      if (nuc->name_ == name){
        return nuc ;
      }
    }
    std::string err_msg = fmt::format("Nuclide {} does not exist in the chain", name) ;
    throw std::invalid_argument(err_msg) ;
}

int Chain::nuclide_index(int nucid) const{
  for (int i = 0; i<nuclides_.size(); i++){
      if (nuclides_[i]->zam_ == nucid) return i ;
  }
  std::string err_msg = fmt::format("Nuclide zam {} does not exist in the chain", nucid) ;
  throw std::invalid_argument(err_msg) ;
}

int Chain::nuclide_index(std::string name) const{
  for (int i = 0; i<nuclides_.size(); i++){
      if (nuclides_[i]->name_ == name) return i ;
  }
  std::string err_msg = fmt::format("Nuclide {} does not exist in the chain", name) ;
  throw std::invalid_argument(err_msg) ;
}

Eigen::SparseMatrix<double> Chain::decayMatrix() const{
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(10000) ;
  for (int inuc = 0; inuc < this->nuclides_.size(); inuc++){
    NuclidePtr nuc = this->nuclides_[inuc] ;
    triplets.push_back(Eigen::Triplet<double>(inuc, inuc, -nuc->dconst_)) ;
    for (const auto& d : nuc->decays_){
      if (!d->targetName_.empty()){
        int jnuc = this->nuclide_index(d->targetName_) ;
        triplets.push_back(Eigen::Triplet<double>(jnuc, inuc, nuc->dconst_ * d->branchingRatio_)) ;
      }
    }
  }
  Eigen::SparseMatrix<double> M(this->nuclides_.size(), this->nuclides_.size()) ;
  M.setFromTriplets(triplets.begin(), triplets.end()) ;
  M.makeCompressed() ;
  return M ;
}

void Chain::removeNuclide(std::string nuc){

}

void Chain::restrict(std::vector<std::string> allowed){
  // std::vector<NuclidePtr> newNuclides ;
  // for (auto& nuc : nuclides_){
  //   if (std::find(allowed.begin(), allowed.end(), nuc->name_) != allowed.end()){
  //     newNuclides.push_back(nuc) ;
  //   }
  // }

  // std::vector<DecayPtr> newDecays ;
  // for (auto& decay : decays_){
  //   bool parentOk = std::find(allowed.begin(), allowed.end(), decay->parentName_) != allowed.end() ;
  //   bool targetOk = std::find(allowed.begin(), allowed.end(), decay->targetName_) != allowed.end() ;
  //   if (parentOk && targetOk){
  //     newDecays.push_back(decay) ;
  //   }
  // }
}