#include "decaysolver.h"
#include <Eigen/Sparse>
#include <cmath>
#include <fmt/format.h>
#include "utils.h"

using msd = std::map<size_t, double> ;

DecaySolver::DecaySolver()= default;

void DecaySolver::compute_coeffs(Chain &chain,
                                 std::map<std::string, double> ccMap) {
  if (!chain.topological_sort())
    throw std::invalid_argument("Chain cannot be topologically sorted");
  chain.tweak_dconst();
  Eigen::SparseMatrix<double> matrix = chain.decayMatrix();

  for (size_t inuc = 0; inuc < chain.nuclides_.size(); inuc++) {
    const auto nuclide = chain.nuclides_[inuc];
    double Cii = nuclide->dconst_;

    if (Cii != 0.) { // UNSTABLE NUCLIDE CASE
      // COMPUTING Ns(i)
      for (const auto &decay : nuclide->decaysUp_) {
        const double br = decay->branchingRatio_;
        const double dconst = decay->parent_->dconst_;
        const double Cik = br * dconst;
        Ns[nuclide->idInChain] += Cik * Ns[decay->parent_->idInChain] / Cii;
      }

      // COMPUTING Fik
      for (size_t knuc = 0; knuc < inuc; knuc++) {
        double Ckk = chain.nuclides_[knuc]->dconst_;
        double factor = 1 / (Cii - Ckk);
        for (const auto &decay : nuclide->decaysUp_) {
          size_t jnuc = decay->parent_->idInChain;
          if (jnuc >= knuc) {
            const double Cij = decay->parent_->dconst_ * decay->branchingRatio_;
            auto jmap = atWithDefault<size_t, msd>(F, jnuc, msd()) ;
            double kval = atWithDefault<size_t, double>(jmap, knuc, 0.) ;
            if (const double val = Cij * kval * factor; val != 0.)
              F[inuc][knuc] += val;
          }
        }
      }

      // COMPUTING Fii
      if (double val = ccMap[nuclide->name_] - Ns[nuclide->idInChain]; val != 0.)
        F[inuc][inuc] = val;
      for (size_t jnuc = 0; jnuc < inuc; jnuc++) {
        auto imap = atWithDefault<size_t, msd>(F, inuc, msd()) ;
        double jval = atWithDefault<size_t, double>(imap, jnuc, 0.) ;
        if (jval != 0.)
          F[inuc][inuc] -= jval;
      }
    } else { // STABLE NUCLIDE CASE
      // COMPUTING Fii
      for (const auto &decay : nuclide->decaysUp_) {
        size_t jnuc = decay->parent_->idInChain;
        double Cij = decay->branchingRatio_ * decay->parent_->dconst_;
        if (double val = Cij * Ns[jnuc]; val != 0.)
          F[inuc][inuc] += val;
      }

      // COMPUTING Fik
      for (size_t knuc = 0; knuc < inuc; knuc++) {
        double Ckk = chain.nuclides_[knuc]->dconst_;
        if (Ckk == 0)
          continue;
        double factor = 1 / (Cii - Ckk);
        for (const auto &decay : nuclide->decaysUp_) {
          size_t jnuc = decay->parent_->idInChain;
          if (jnuc >= knuc) {
            double Cij = decay->parent_->dconst_ * decay->branchingRatio_;
            auto jmap = atWithDefault<size_t, msd>(F, jnuc, msd()) ;
            double kval = atWithDefault<size_t, double>(jmap, knuc, 0.) ;
            if (kval != 0.)
              F[inuc][knuc] += Cij * kval * factor;
          }
        }
      }

      // COMPUTING Ns(i)
      Ns[inuc] = ccMap[nuclide->name_];
      for (size_t knuc = 0; knuc < inuc; knuc++) {
        auto imap = atWithDefault<size_t, msd>(F, inuc, msd()) ;
        double kval = atWithDefault<size_t, double>(imap, knuc, 0.) ;
        if (kval != 0.)
          Ns[inuc] -= kval;
      }
    }
  }
}

std::vector<std::vector<double>>
DecaySolver::run(Chain &chain, const std::map<std::string, double>& ccMap,
                 std::vector<double> times) {
  const size_t nt = times.size();
  const size_t nn = chain.nuclides_.size();
  this->compute_coeffs(chain, ccMap);
  std::vector<std::vector<double>> N(nt, std::vector<double>(nn, 0));
  for (auto const &[key, val] : ccMap) {
    const size_t inuc = chain.nuclide_index(key);
    N[0][inuc] = val;
  }

  for (size_t it = 1; it < nt; it++) {
    fmt::print("{:.4e} -> {:.4e}\n", times[it-1], times[it]);
    for (size_t inuc = 0; inuc < nn; inuc++) {
      auto imap = atWithDefault(F, inuc, msd());
      N[it][inuc] += Ns[inuc];
      for (const auto& [jnuc, Fij]: imap) {
        if (inuc == jnuc && chain.nuclides_[inuc]->dconst_ == 0)
          continue;
        N[it][inuc] += Fij * std::exp(-chain.nuclides_[jnuc]->dconst_ * times[it]);
      }
      if (chain.nuclides_[inuc]->dconst_ == 0)
        N[it][inuc] += atWithDefault(imap, inuc, 0.) * times[it] ;
    }
  }

  std::vector<std::string> nuclidenames;
  for (const auto &nuclide : chain.nuclides_)
    nuclidenames.push_back(nuclide->name_);
  fmt::print("times.size()={}\n", times.size()) ;
  fmt::print("nuclidenames.size()={}\n", nuclidenames.size()) ;
  fmt::print("N.size()={}\n", N.size()) ;
  fmt::print("N[0].size()={}\n", N[0].size()) ;
  results_ = Results(N, nuclidenames, times);
  return N;
}
