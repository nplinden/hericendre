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

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

enum class RunMode
{
  PROCESS,
  DEPLETE,
  BURN
};

int deplete()
{
  Chain chain("/home/nlinden/workspace/deplete/data/chain-jeff33-long.xml");
  std::map<std::string, double> ccMap;

  Eigen::VectorXd N(chain.nuclides_.size());
  for (int i = 0; i < chain.nuclides_.size(); i++)
    N(i) = 1.;
  Solver solver;
  Eigen::VectorXd res = solver.run(chain, N, 760854000000);
  return 0;
}

int display_info()
{
  Chain chain("/home/nlinden/workspace/deplete/data/chain-jeff33-long.xml");
  // Chain chain("/home/nlinden/workspace/deplete/data/short.xml") ;
  std::map<std::string, double> ccMap;

  fmt::print("Number of nuclides: {}\n", chain.nuclides_.size());
  fmt::print("Number of decays: {}\n", chain.decays_.size());
  fmt::print("Number of nreactions: {}\n", chain.reactions_.size());

  auto nuc = chain.find(922350);
  fmt::print("{}\n", *nuc);
  return 1;
}

void process(int argc, char *argv[])
{
  char *inFile;
  char *outFile;
  for (int count = 0; count < argc; ++count)
  {
    if (std::string(argv[count]) == "-c")
    {
      inFile = argv[count + 1];
    }
    if (std::string(argv[count]) == "-o")
    {
      outFile = argv[count + 1];
    }
  }
  Chain chain(inFile);
  chain.write(outFile);
}

void info(int argc, char *argv[])
{
  char *inFile;
  for (int count = 0; count < argc; ++count)
  {
    if (std::string(argv[count]) == "-c")
    {
      inFile = argv[count + 1];
    }
  }
}

void dfs(int nucid, std::vector<bool> &visited, Chain &chain)
{
  if (visited[nucid])
    return;
  visited[nucid] = true;

  std::vector<DecayPtr> decays = chain.nuclides_[nucid]->decays_;
  std::vector<int> neighbours;
  fmt::print("{}: {}\n", nucid, chain.nuclides_[nucid]->name_);
  for (DecayPtr decay : chain.nuclides_[nucid]->decays_)
  {
    fmt::print("\t{}\n", decay->target_->name_);
    neighbours.push_back(decay->target_->idInChain);
  }
  for (const auto neighboursId : neighbours)
  {
    dfs(neighboursId, visited, chain);
  }
}

void reachable(Chain chain, std::string nucname)
{
  int n = chain.nuclides_.size();
  std::vector<bool> visited;
  for (int i = 0; i < n; i++)
    visited.push_back(false);
  int initialId = chain.find(nucname)->idInChain;
  dfs(initialId, visited, chain);

  for (int i = 0; i < n; i++)
  {
    // if (visited[i])
    // fmt::print("{}\n", chain.nuclides_[i]->name_) ;
  }
}

int main(int argc, char *argv[])
{

  Chain chain("/home/nlinden/workspace/hericendre/data/chain_endfb71_sfr.xml");
  reachable(chain, "Pu239");

  // RunMode mode ;

  if (argc == 1)
  {
    fmt::print("Please provide instructions\n");
    return 0;
  }
  else
  {
    if (std::string(argv[1]) == "process")
    {
      process(argc, argv);
    }
    if (std::string(argv[1]) == "info")
    {
      info(argc, argv);
    }
    return 0;
  }

  return 0;
}
