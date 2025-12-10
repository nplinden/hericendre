#include "decaysolver.h"
#include <Eigen/Sparse>
#include <cmath>
#include <fmt/format.h>
#include "utils.h"

using msd = std::map<size_t, double>;

DecaySolver::DecaySolver(Chain &chain)
{
  chain_ = chain;
};

void DecaySolver::compute_coeffs(std::map<std::string, double> ccMap)
{
  if (!chain_.topological_sort())
  {
    throw std::invalid_argument("Chain cannot be topologically sorted");
  }
  chain_.tweak_dconst();
  Eigen::SparseMatrix<double> matrix = chain_.decayMatrix();

  for (const NuclidePtr nuclide : chain_.nuclides_)
  {
    compute_Fik(nuclide);
    compute_Ns(nuclide, ccMap);
    compute_Fii(nuclide, ccMap);
  }
}

Results
DecaySolver::run(const std::map<std::string, double> &ccMap,
                 std::vector<double> times)
{
  const size_t nt = times.size();
  const size_t nn = chain_.nuclides_.size();
  this->compute_coeffs(ccMap);
  std::vector<std::vector<double>> N(nt, std::vector<double>(nn, 0));
  for (auto const &[key, val] : ccMap)
  {
    const size_t i = chain_.nuclide_index(key);
    N[0][i] = val;
  }

  for (size_t it = 1; it < nt; it++)
  {
    fmt::print("{:.4e} -> {:.4e}\n", times[it - 1], times[it]);
    for (size_t i = 0; i < nn; i++)
    {
      N[it][i] += Ns[i];

      auto Fi = F.find(i);
      if (Fi == F.end())
        continue;

      for (const auto &[j, Fij] : Fi->second)
      {
        if (i == j && chain_.nuclides_[i]->dconst_ == 0)
          continue;
        N[it][i] += Fij * std::exp(-chain_.nuclides_[j]->dconst_ * times[it]);
      }

      auto Fii_it = Fi->second.find(i);
      if (Fii_it != Fi->second.end())
        continue;

      double Fii = Fii_it->second;
      if (chain_.nuclides_[i]->dconst_ == 0)
        N[it][i] += Fii * times[it];
    }
  }

  std::vector<std::string> nuclidenames;
  for (const auto &nuclide : chain_.nuclides_)
    nuclidenames.push_back(nuclide->name_);
  return Results(N, nuclidenames, times);
}

void DecaySolver::compute_Fik(const NuclidePtr nuclide)
{
  const size_t i = nuclide->idInChain;
  const double Cii = nuclide->dconst_;

  for (size_t k = 0; k < i; k++)
  {
    double Ckk = chain_.nuclides_[k]->dconst_;

    if (Ckk == 0 && Cii == 0)
      continue;

    double factor = 1 / (Cii - Ckk); // Cii != Ckk since tweak_dconst was called
    for (const auto &decay : nuclide->decaysUp_)
    {
      size_t j = decay->parent_->idInChain;
      if (j >= k)
      {
        const double Cij = decay->parent_->dconst_ * decay->branchingRatio_;

        auto Fj = F.find(j);
        if (Fj == F.end())
          continue;

        auto Fjk_it = Fj->second.find(k);
        if (Fjk_it == Fj->second.end())
          continue;

        double Fjk = Fjk_it->second;

        double val = Cij * Fjk * factor;
        if (val != 0.)
          F[i][k] += val;
      }
    }
  }
}

void DecaySolver::compute_Ns(const NuclidePtr nuclide, const std::map<std::string, double> &ccMap)
{
  if (!nuclide->isStable())
  {
    /*Ns = 0 for unstable nuclides, when there are no external sources
    TODO: add external sources capability
     */
    return;
  }
  const size_t i = nuclide->idInChain;
  double N0 = (ccMap.find(nuclide->name_) != ccMap.end()) ? ccMap.at(nuclide->name_) : 0.0;

  Ns[i] = N0;

  auto Fi = F.find(i);
  if (Fi == F.end())
    return;

  for (const auto &[k, Fik] : Fi->second)
  {
    if (k == i)
      continue;
    if (Fik == 0.)
      continue;
    Ns[i] -= Fik;
  }
}

void DecaySolver::compute_Fii(const NuclidePtr nuclide, const std::map<std::string, double> &ccMap)
{
  if (nuclide->isStable())
  {
    /* For stable nuclides, Fii is defined as:
      Fii = sum_{j} lambda(j) * g(j -> i) * Ns(j)

      In this sum, either:
        - j is unstable, in which case Ns(j) = 0
        - j is stable, in which case lambda(j) = 0
      Therefore, Fii = 0 for stable nuclides.

    Note that Ns(j) = 0 for unstable nuclides because there are no external sources
    TODO: add external sources capability
    */
    return;
  }
  const size_t i = nuclide->idInChain;

  double N0 = (ccMap.find(nuclide->name_) != ccMap.end()) ? ccMap.at(nuclide->name_) : 0.0;

  if (double val = N0 - Ns[nuclide->idInChain]; val != 0.)
    F[i][i] = val;
  for (size_t j = 0; j < i; j++)
  {
    auto Fi = F.find(i);
    if (Fi == F.end())
      continue;

    auto Fij_it = Fi->second.find(j);
    if (Fij_it == Fi->second.end())
      continue;

    double Fij = Fij_it->second;
    if (Fij != 0.)
      F[i][i] -= Fij;
  }
}