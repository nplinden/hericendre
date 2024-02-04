#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "pugixml.hpp"
#include <fmt/core.h>
#include <fmt/os.h>
#include <nuclide.hpp>
#include <decay.hpp>
#include <chain.hpp>
#include <Eigen/Sparse>
#include <map>
#include <solver.hpp>
#include <fmt/ranges.h>

typedef Eigen::SparseMatrix<double> SpMat ;
typedef Eigen::Triplet<double> T;

int deplete(){
  Chain chain("/home/nlinden/workspace/deplete/data/chain-jeff33-long.xml") ;
  std::map<std::string, double> ccMap ;

  Eigen::VectorXd N(chain.nuclides_.size()) ;
  for (int i = 0; i < chain.nuclides_.size(); i++)
    N(i) = 1. ;
  Solver solver;
  Eigen::VectorXd res = solver.run(chain, N, 760854000000) ;
  return 0 ;
}

int display_info(){
  Chain chain("/home/nlinden/workspace/deplete/data/chain-jeff33-long.xml") ;
  // Chain chain("/home/nlinden/workspace/deplete/data/short.xml") ;
  std::map<std::string, double> ccMap ;


  fmt::print("Number of nuclides: {}\n", chain.nuclides_.size()) ;
  fmt::print("Number of decays: {}\n", chain.decays_.size()) ;
  fmt::print("Number of nreactions: {}\n", chain.nreactions_.size()) ;

  auto nuc = chain.find(922350) ;
  fmt::print("{}\n", *nuc) ;
  return 1 ;
}

int main(int argc, char* argv[])
{
  Chain chain("/home/nlinden/workspace/deplete/data/chain-jeff33-long.xml") ;
  // Chain chain("/home/nlinden/workspace/deplete/data/short.xml") ;

  pugi::xml_document doc;
  auto root = doc.append_child("depletion_chain");

  for (auto nuclide : chain.nuclides_) {
    nuclide->addNode(root) ;
  }

  doc.save_file("chain.xml", "  ") ;

  // auto nuc = chain.find("Pu239") ;
  // for (auto s : nuc->sources_){
  //   fmt::print("{} {}\n", s->particle_, s->type_) ;
  //   if (s->type_ != "mixture"){
  //     fmt::print("{}\n", s->energy_.size()) ;
  //     fmt::print("{}\n", s->intensity_.size()) ;
  //   }
  // }

  return 0;
}
