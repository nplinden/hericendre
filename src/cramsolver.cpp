#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <fmt/core.h>
#include <fmt/os.h>
#include <cramsolver.h>

CRAMSolver::CRAMSolver() = default;

std::vector<Eigen::VectorXd> CRAMSolver::run(const Chain &chain,
                                             const Eigen::VectorXd &ccVector, const double dt) const {
    const size_t n_nuclides = chain.nuclides_.size();
    const SpComplex M = chain.decayMatrix().cast<cdouble>();

    std::vector<TrComplex> t_triplets;
    t_triplets.reserve(n_nuclides);
    for (size_t j = 0; j < n_nuclides; j++)
        t_triplets.emplace_back(j, j, std::complex<double>(1., 0.));
    SpComplex Identity(n_nuclides, n_nuclides);
    Identity.setFromTriplets(t_triplets.begin(), t_triplets.end());

    Eigen::VectorX<cdouble> N = ccVector.cast<cdouble>();
    Eigen::SparseLU<SpComplex, Eigen::COLAMDOrdering<int> > solver;
    for (size_t i = 0; i < theta48.size(); i++) {
        cdouble theta = theta48[i];
        cdouble alpha = alpha48[i];

        SpComplex A = M * dt - theta * Identity;
        solver.analyzePattern(A);
        solver.factorize(A);
        Eigen::VectorX<cdouble> x = (alpha * solver.solve(N));
        N += 2 * x.real();
    }

    Eigen::VectorX<double> realN = alpha48_0 * N.real();

    for (double &i: realN) {
        if (i < cutoff_)
            i = 0.;
    }

    return std::vector<Eigen::VectorXd>({ccVector, realN});
}

std::vector<Eigen::VectorXd> CRAMSolver::run(const Chain &chain,
                                             const Eigen::VectorXd &ccVector,
                                             std::vector<double> times) {
    std::vector<Eigen::VectorXd> concentrations;
    Eigen::VectorXd N(ccVector);

    concentrations.push_back(N);
    for (size_t it = 1; it < times.size(); it++) {
        const double dt = times[it] - times[it - 1];
        fmt::print("{:.4e} -> {:.4e}\n", times[it - 1], times[it]);
        N = run(chain, N, dt).back();
        concentrations.push_back(N);
    }

    std::vector<std::string> nuclides;
    for (const auto &nuclide: chain.nuclides_) {
        nuclides.push_back(nuclide->name_);
    }
    results_ = Results(concentrations, nuclides, times);
    return concentrations;
}

std::vector<Eigen::VectorXd> CRAMSolver::run(const Chain &chain,
                                             const std::map<std::string, double> &ccMap,
                                             const std::vector<double> &times) {
    Eigen::VectorXd N(chain.nuclides_.size());
    for (size_t i = 0; i < chain.nuclides_.size(); i++)
        N(i) = 0.;

    for (auto const &[key, val]: ccMap) {
        const size_t inuc = chain.nuclide_index(key);
        N(inuc) = val;
    }

    return run(chain, N, times);
}
