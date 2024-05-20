#include "decaysolver.h"
#include <Eigen/Sparse>
#include <cmath>

DecaySolver::DecaySolver(){};

void DecaySolver::compute_coeffs(Chain &chain,
                                 std::map<std::string, double> ccMap) {
  bool success = chain.topological_sort();
  if (!success)
    throw std::invalid_argument("Chain cannot be topologically sorted");
  chain.tweak_dconst();
  Eigen::SparseMatrix<double> matrix = chain.decayMatrix();

  for (size_t inuc = 0; inuc < chain.nuclides_.size(); inuc++) {
    auto nuclide = chain.nuclides_[inuc];
    // fmt::print("i={} {}\n", inuc, nuclide->name_);
    // fmt::print("\tλ={}\n", nuclide->dconst_);
    double Cii = nuclide->dconst_;
    // fmt::print(" Cii={}\n", Cii);

    // UNSTABLE NUCLIDE CASE
    if (Cii != 0.) {
      // COMPUTING Ns(i)
      for (const auto &decay : nuclide->decaysUp_) {
        double br = decay->branchingRatio_;
        double dconst = decay->parent_->dconst_;
        double Cik = br * dconst;
        Ns[nuclide->idInChain] += Cik * Ns[decay->parent_->idInChain] / Cii;
      }
      // fmt::print("\tNs({})={}\n", nuclide->name_, Ns[nuclide->idInChain]);

      // COMPUTING Fik
      for (size_t knuc = 0; knuc < inuc; knuc++) {
        double Ckk = chain.nuclides_[knuc]->dconst_;
        double factor = 1 / (Cii - Ckk);
        for (const auto &decay : nuclide->decaysUp_) {
          size_t jnuc = decay->parent_->idInChain;
          if (jnuc >= knuc) {
            double Cij = decay->parent_->dconst_ * decay->branchingRatio_;
            F[inuc][knuc] += Cij * F[jnuc][knuc] * factor;
          }
        }
        // if (F[inuc][knuc] != 0)
        // fmt::print("\tF[{}][{}]={}\n", nuclide->name_,
        // chain.nuclides_[knuc]->name_, F[inuc][knuc]);
      }

      // COMPUTING Fii
      F[inuc][inuc] = ccMap[nuclide->name_] - Ns[nuclide->idInChain];
      for (size_t jnuc = 0; jnuc < inuc; jnuc++) {
        F[inuc][inuc] -= F[inuc][jnuc];
      }
      // fmt::print("\tF[{}][{}]={}\n", nuclide->name_, nuclide->name_,
      // F[inuc][inuc]);
      // STABLE NUCLIDE CASE
    } else {
      // COMPUTING Fii
      for (const auto &decay : nuclide->decaysUp_) {
        size_t jnuc = decay->parent_->idInChain;
        double Cij = decay->branchingRatio_ * decay->parent_->dconst_;
        F[inuc][inuc] += Cij * Ns[jnuc];
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
            F[inuc][knuc] += Cij * F[jnuc][knuc] * factor;
          }
        }
        // if (F[inuc][knuc] != 0)
        // fmt::print("\tF[{}][{}]={}\n", nuclide->name_,
        //            chain.nuclides_[knuc]->name_, F[inuc][knuc]);
      }

      // COMPUTING Ns(i)
      Ns[inuc] = ccMap[nuclide->name_];
      for (size_t knuc = 0; knuc < inuc; knuc++) {
        Ns[inuc] -= F[inuc][knuc];
      }
      // fmt::print("\tNs({})={}\n", nuclide->name_, Ns[nuclide->idInChain]);
    }
  }
}

std::vector<std::vector<double>>
DecaySolver::run(Chain &chain, std::map<std::string, double> ccMap,
                 std::vector<double> times) {
  size_t nt = times.size();
  size_t nn = chain.nuclides_.size();
  this->compute_coeffs(chain, ccMap);
  std::vector<std::vector<double>> N(nt + 1, std::vector<double>(nn, 0));
  for (auto const &[key, val] : ccMap) {
    int inuc = chain.nuclide_index(key);
    N[0][inuc] = val;
  }

  for (size_t it = 0; it < nt; it++) {
    fmt::print("{} {}\n", it, times[it]);
    for (size_t inuc = 0; inuc < nn; inuc++) {
      auto nuclide = chain.nuclides_[inuc];
      if (nuclide->dconst_ != 0) {
        // unstable case
        N[it + 1][inuc] += Ns[inuc];
        for (size_t jnuc = 0; jnuc <= inuc; jnuc++)
          N[it + 1][inuc] +=
              F[inuc][jnuc] *
              std::exp(-chain.nuclides_[jnuc]->dconst_ * times[it]);
      } else { // stable case
        N[it + 1][inuc] += Ns[inuc];
        for (size_t jnuc = 0; jnuc < inuc; jnuc++)
          N[it + 1][inuc] +=
              F[inuc][jnuc] *
              std::exp(-chain.nuclides_[jnuc]->dconst_ * times[it]);
        N[it + 1][inuc] += F[inuc][inuc] * times[it];
      }
    }
  }

  std::vector<double> times_with_zero;
  times_with_zero.push_back(0.);
  for (const auto &t : times) {
    times_with_zero.push_back(t);
  }
  std::vector<std::string> nuclidenames;
  for (const auto &nuclide : chain.nuclides_)
    nuclidenames.push_back(nuclide->name_);
  results_ = Results(N, nuclidenames, times_with_zero);
  return N;
}
