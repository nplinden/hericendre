#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <fmt/core.h>
#include <fmt/os.h>
#include <cramsolver.h>

CRAMSolver::CRAMSolver() = default;

Eigen::VectorXd CRAMSolver::run(const SpComplex &M, const Eigen::VectorXd &ccVector, const double dt) const {
    const size_t n = M.rows();
    // const SpComplex M = chain.decayMatrix().cast<cdouble>();

    SpComplex Identity(n, n);
    Identity.setIdentity();

    Eigen::VectorX<cdouble> N = ccVector.cast<cdouble>();
    Eigen::SparseLU<SpComplex, Eigen::COLAMDOrdering<int> > solver;

    SpComplex A_pattern = M * dt - theta48[0] * Identity;
    solver.analyzePattern(A_pattern);

    for (size_t i = 0; i < theta48.size(); i++) {
        const cdouble theta = theta48[i];
        const cdouble alpha = alpha48[i];

        SpComplex A = M * dt - theta * Identity;
        solver.factorize(A);
        Eigen::VectorX<cdouble> x = alpha * solver.solve(N);
        N += 2 * x.real();
    }

    Eigen::VectorX<double> realN = alpha48_0 * N.real();

    for (double &val: realN) {
        if (val < cutoff_) {
            val = 0.;
        }
    }

    return realN;
}

Results CRAMSolver::run(const Chain &chain,
                                             const Eigen::VectorXd &ccVector,
                                             const std::vector<double> &times) {
    const SpComplex M = chain.decayMatrix().cast<cdouble>();

    std::vector<Eigen::VectorXd> concentrations;
    concentrations.reserve(times.size());

    Eigen::VectorXd N = ccVector;
    concentrations.push_back(N);

    for (size_t it = 1; it < times.size(); it++) {
        const double dt = times[it] - times[it - 1];
        fmt::print("{:.4e} -> {:.4e}\n", times[it - 1], times[it]);
        N = run(M, N, dt);
        concentrations.push_back(N);
    }

    return Results(concentrations, chain.name_vector(), times);
}

Results CRAMSolver::run(const Chain &chain,
                                             const std::map<std::string, double> &ccMap,
                                             const std::vector<double> &times) {
    Eigen::VectorXd N = Eigen::VectorXd::Zero(chain.nuclides_.size());

    for (const auto &[key, val]: ccMap) {
        const size_t inuc = chain.nuclide_index(key);
        N(inuc) = val;
    }

    return run(chain, N, times);
}
