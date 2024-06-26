#include "decaysolver.h"
#include <Eigen/Sparse>
#include <cmath>
#include <fmt/format.h>
#include "utils.h"

using msd = std::map<size_t, double>;

DecaySolver::DecaySolver(Chain& chain) {
    chain_ = chain;
};

void DecaySolver::compute_coeffs(std::map<std::string, double> ccMap) {
  ccMap_ = ccMap;
  if (!chain_.topological_sort())
    throw std::invalid_argument("Chain cannot be topologically sorted");
  chain_.tweak_dconst();
  Eigen::SparseMatrix<double> matrix = chain_.decayMatrix();

  for (size_t inuc = 0; inuc < chain_.nuclides_.size(); inuc++) {
    const auto nuclide = chain_.nuclides_[inuc];
    double Cii = nuclide->dconst_;

    if (Cii != 0.) {
      // UNSTABLE NUCLIDE CASE
      // COMPUTING Ns(i)
      for (const auto &decay: nuclide->decaysUp_) {
        const double br = decay->branchingRatio_;
        const double dconst = decay->parent_->dconst_;
        const double Cik = br * dconst;
        Ns[nuclide->idInChain] += Cik * Ns[decay->parent_->idInChain] / Cii;
      }

      // COMPUTING Fik
      for (size_t knuc = 0; knuc < inuc; knuc++) {
        double Ckk = chain_.nuclides_[knuc]->dconst_;
        double factor = 1 / (Cii - Ckk);
        for (const auto &decay: nuclide->decaysUp_) {
          size_t jnuc = decay->parent_->idInChain;
          if (jnuc >= knuc) {
            const double Cij = decay->parent_->dconst_ * decay->branchingRatio_;
            auto jmap = atWithDefault<size_t, msd>(F, jnuc, msd());
            double kval = atWithDefault<size_t, double>(jmap, knuc, 0.);
            if (const double val = Cij * kval * factor; val != 0.)
              F[inuc][knuc] += val;
          }
        }
      }

      // COMPUTING Fii
      if (double val = ccMap[nuclide->name_] - Ns[nuclide->idInChain]; val != 0.)
        F[inuc][inuc] = val;
      for (size_t jnuc = 0; jnuc < inuc; jnuc++) {
        auto imap = atWithDefault<size_t, msd>(F, inuc, msd());
        double jval = atWithDefault<size_t, double>(imap, jnuc, 0.);
        if (jval != 0.)
          F[inuc][inuc] -= jval;
      }
    } else {
      // STABLE NUCLIDE CASE
      // COMPUTING Fii
      for (const auto &decay: nuclide->decaysUp_) {
        size_t jnuc = decay->parent_->idInChain;
        double Cij = decay->branchingRatio_ * decay->parent_->dconst_;
        if (double val = Cij * Ns[jnuc]; val != 0.)
          F[inuc][inuc] += val;
      }

      // COMPUTING Fik
      for (size_t knuc = 0; knuc < inuc; knuc++) {
        double Ckk = chain_.nuclides_[knuc]->dconst_;
        if (Ckk == 0)
          continue;
        double factor = 1 / (Cii - Ckk);
        for (const auto &decay: nuclide->decaysUp_) {
          size_t jnuc = decay->parent_->idInChain;
          if (jnuc >= knuc) {
            double Cij = decay->parent_->dconst_ * decay->branchingRatio_;
            auto jmap = atWithDefault<size_t, msd>(F, jnuc, msd());
            double kval = atWithDefault<size_t, double>(jmap, knuc, 0.);
            if (kval != 0.)
              F[inuc][knuc] += Cij * kval * factor;
          }
        }
      }

      // COMPUTING Ns(i)
      Ns[inuc] = ccMap[nuclide->name_];
      for (size_t knuc = 0; knuc < inuc; knuc++) {
        auto imap = atWithDefault<size_t, msd>(F, inuc, msd());
        double kval = atWithDefault<size_t, double>(imap, knuc, 0.);
        if (kval != 0.)
          Ns[inuc] -= kval;
      }
    }
  }
}

std::vector<std::vector<double> >
DecaySolver::run(const std::map<std::string, double> &ccMap,
                 std::vector<double> times) {
  const size_t nt = times.size();
  const size_t nn = chain_.nuclides_.size();
  this->compute_coeffs(ccMap);
  std::vector<std::vector<double> > N(nt, std::vector<double>(nn, 0));
  for (auto const &[key, val]: ccMap) {
    const size_t inuc = chain_.nuclide_index(key);
    N[0][inuc] = val;
  }

  for (size_t it = 1; it < nt; it++) {
    fmt::print("{:.4e} -> {:.4e}\n", times[it - 1], times[it]);
    for (size_t inuc = 0; inuc < nn; inuc++) {
      auto imap = atWithDefault(F, inuc, msd());
      N[it][inuc] += Ns[inuc];
      for (const auto &[jnuc, Fij]: imap) {
        if (inuc == jnuc && chain_.nuclides_[inuc]->dconst_ == 0)
          continue;
        N[it][inuc] += Fij * std::exp(-chain_.nuclides_[jnuc]->dconst_ * times[it]);
      }
      if (chain_.nuclides_[inuc]->dconst_ == 0)
        N[it][inuc] += atWithDefault(imap, inuc, 0.) * times[it];
    }
  }

  std::vector<std::string> nuclidenames;
  for (const auto &nuclide: chain_.nuclides_)
    nuclidenames.push_back(nuclide->name_);
  fmt::print("times.size()={}\n", times.size());
  fmt::print("nuclidenames.size()={}\n", nuclidenames.size());
  fmt::print("N.size()={}\n", N.size());
  fmt::print("N[0].size()={}\n", N[0].size());
  results_ = Results(N, nuclidenames, times);
  return N;
}

void DecaySolver::to_hdf5(H5Easy::File &file) const {
  for (auto const& [i, dico]: this->F){
      std::string inuc = this->chain_.nuclides_[i]->name_;
      for (auto const& [j, val]: dico){
          std::string jnuc = this->chain_.nuclides_[j]->name_;
          std::string h5path = fmt::format("/solver/F/{}/{}", inuc, jnuc);
          H5Easy::dump(file, h5path, std::vector<double>({val}));
      }
  }
  H5Easy::dump(file, "/solver/dconst", this->chain_.dconst_vector());
  H5Easy::dump(file, "/solver/nuclides", this->chain_.name_vector());
  H5Easy::dump(file, "/solver/Ns", this->Ns_vector());
}

std::vector<double> DecaySolver::Ns_vector() const {
    std::vector<double> vec;
    for (size_t i = 0; i < this->chain_.nuclides_.size(); i++){
       vec.push_back(this->Ns.at(i));
    }
    return vec;
}

