#ifndef SOLVER_HPP_INCLUDED
#define SOLVER_HPP_INCLUDED
#include <Eigen/Sparse>
#include <chain.h>
#include <microxs.h>
#include <complex>
#include <map>
#include <results.h>
#include <vector>

using cdouble = std::complex<double>;
using SpComplex = Eigen::SparseMatrix<cdouble>;
using TrComplex = Eigen::Triplet<cdouble>;
using MicroXSPtr = std::shared_ptr<MicroXS>;

/**
 * @brief Chebyshev Rational Approximation Method (CRAM) solver for depletion equations
 *
 * This class implements the CRAM-48 algorithm for solving the Bateman equations,
 * which describe radioactive decay and transmutation chains. CRAM-48 uses a 48th-order
 * rational approximation with complex poles to accurately compute the matrix exponential
 * exp(M*t), where M is the decay/transmutation matrix.
 *
 * The CRAM method is highly accurate for stiff systems and is particularly well-suited
 * for nuclear depletion problems where decay constants can span many orders of magnitude.
 *
 * Reference: M. Pusa, "Rational Approximations to the Matrix Exponential in Burnup Calculations"
 */
class CRAMSolver
{
public:
    CRAMSolver(MicroXSPtr microxs): microxs_(microxs) {};

    /**
     * @brief Solve depletion for multiple time steps with initial concentrations from a map
     *
     * Solves the Bateman equations from t=0 through all provided time points.
     * Initial concentrations are specified by nuclide name. Nuclides not in the map
     * are initialized to zero.
     *
     * @param chain The depletion chain containing nuclides and decay data
     * @param ccMap Initial concentrations as a map: ccMap["U235"] = 1.0
     * @param times Time points (in seconds) at which to compute concentrations
     * @return Results object containing concentrations at all time points
     */
    Results run(const Chain &chain,
                const std::map<std::string, double> &ccMap,
                const std::vector<double> &times);

    /**
     * @brief Solve depletion for multiple time steps with initial concentrations from a vector
     *
     * Solves the Bateman equations from t=0 through all provided time points.
     * Initial concentrations are provided as an Eigen vector matching the chain nuclide order.
     *
     * @param chain The depletion chain containing nuclides and decay data
     * @param ccVector Initial concentration vector (length must match chain size)
     * @param times Time points (in seconds) at which to compute concentrations
     * @return Results object containing concentrations at all time points
     */
    Results run(const Chain &chain, const Eigen::VectorXd &ccVector,
                const std::vector<double> &times);

    /**
     * @brief Concentration cutoff threshold for numerical stability
     *
     * Concentrations below this value are set to zero after each time step.
     * This prevents numerical noise from accumulating in very small concentrations.
     * Default: 1e-14
     */
    double cutoff_ = 1.e-14;

    std::shared_ptr<MicroXS> microxs_;

private:
    /**
     * @brief Core CRAM-48 solver for a single time step
     *
     * Computes the solution to dN/dt = M*N over a single time interval dt
     * using the 48th-order Chebyshev Rational Approximation Method.
     *
     * @param M Pre-computed decay/transmutation matrix (complex sparse)
     * @param ccVector Initial concentration vector at time t
     * @param dt Time step size in seconds
     * @return Concentration vector at time t+dt
     */
    Eigen::VectorXd run(const SpComplex &M, const Eigen::VectorXd &ccVector,
                        double dt) const;

    /**
     * @brief Complex poles (theta_i) for CRAM-48 approximation
     *
     * These are the 24 complex conjugate pole pairs used in the partial fraction
     * expansion of the (16,16) Padé approximation to exp(z).
     */
    std::vector<cdouble> theta48 = {
        cdouble(-4.465731934165702e+1, +6.233225190695437e+1),
        cdouble(-5.284616241568964e+0, +4.057499381311059e+1),
        cdouble(-8.867715667624458e+0, +4.325515754166724e+1),
        cdouble(+3.493013124279215e+0, +3.281615453173585e+1),
        cdouble(+1.564102508858634e+1, +1.558061616372237e+1),
        cdouble(+1.742097597385893e+1, +1.076629305714420e+1),
        cdouble(-2.834466755180654e+1, +5.492841024648724e+1),
        cdouble(+1.661569367939544e+1, +1.316994930024688e+1),
        cdouble(+8.011836167974721e+0, +2.780232111309410e+1),
        cdouble(-2.056267541998229e+0, +3.794824788914354e+1),
        cdouble(+1.449208170441839e+1, +1.799988210051809e+1),
        cdouble(+1.853807176907916e+1, +5.974332563100539e+0),
        cdouble(+9.932562704505182e+0, +2.532823409972962e+1),
        cdouble(-2.244223871767187e+1, +5.179633600312162e+1),
        cdouble(+8.590014121680897e-1, +3.536456194294350e+1),
        cdouble(-1.286192925744479e+1, +4.600304902833652e+1),
        cdouble(+1.164596909542055e+1, +2.287153304140217e+1),
        cdouble(+1.806076684783089e+1, +8.368200580099821e+0),
        cdouble(+5.870672154659249e+0, +3.029700159040121e+1),
        cdouble(-3.542938819659747e+1, +5.834381701800013e+1),
        cdouble(+1.901323489060250e+1, +1.194282058271408e+0),
        cdouble(+1.885508331552577e+1, +3.583428564427879e+0),
        cdouble(-1.734689708174982e+1, +4.883941101108207e+1),
        cdouble(+1.316284237125190e+1, +2.042951874827759e+1)};

    /**
     * @brief Complex residues (alpha_i) for CRAM-48 approximation
     *
     * These are the 24 complex conjugate residue pairs corresponding to the poles
     * in theta48. Together with theta48, they define the rational approximation.
     */
    std::vector<cdouble> alpha48 = {
        cdouble(+6.387380733878774e+2, -6.743912502859256e+2),
        cdouble(+1.909896179065730e+2, -3.973203432721332e+2),
        cdouble(+4.236195226571914e+2, -2.041233768918671e+3),
        cdouble(+4.645770595258726e+2, -1.652917287299683e+3),
        cdouble(+7.765163276752433e+2, -1.783617639907328e+4),
        cdouble(+1.907115136768522e+3, -5.887068595142284e+4),
        cdouble(+2.909892685603256e+3, -9.953255345514560e+3),
        cdouble(+1.944772206620450e+2, -1.427131226068449e+3),
        cdouble(+1.382799786972332e+5, -3.256885197214938e+6),
        cdouble(+5.628442079602433e+3, -2.924284515884309e+4),
        cdouble(+2.151681283794220e+2, -1.121774011188224e+3),
        cdouble(+1.324720240514420e+3, -6.370088443140973e+4),
        cdouble(+1.617548476343347e+4, -1.008798413156542e+6),
        cdouble(+1.112729040439685e+2, -8.837109731680418e+1),
        cdouble(+1.074624783191125e+2, -1.457246116408180e+2),
        cdouble(+8.835727765158191e+1, -6.388286188419360e+1),
        cdouble(+9.354078136054179e+1, -2.195424319460237e+2),
        cdouble(+9.418142823531573e+1, -6.719055740098035e+2),
        cdouble(+1.040012390717851e+2, -1.693747595553868e+2),
        cdouble(+6.861882624343235e+1, -1.177598523430493e+1),
        cdouble(+8.766654491283722e+1, -4.596464999363902e+3),
        cdouble(+1.056007619389650e+2, -1.738294585524067e+3),
        cdouble(+7.738987569039419e+1, -4.311715386228984e+1),
        cdouble(+1.041366366475571e+2, -2.777743732451969e+2)};

    /**
     * @brief Scalar coefficient (alpha_0) for CRAM-48 approximation
     *
     * This is the real scalar term in the CRAM-48 partial fraction expansion.
     * It scales the final result after all complex pole contributions are summed.
     */
    double alpha48_0 = 2.258038182743983e-47;
};

#endif
