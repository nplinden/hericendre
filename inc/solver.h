#ifndef SOLVER_HPP_INCLUDED
#define SOLVER_HPP_INCLUDED
#include <Eigen/Sparse>
#include <chain.h>
#include <complex>
#include <map>
#include <results.h>
#include <vector>

using cdouble = std::complex<double>;
using SpComplex = Eigen::SparseMatrix<cdouble>;
using TrComplex = Eigen::Triplet<cdouble>;

class Solver {
public:
 Solver();

 /**
  * \brief Solve the Bateman equation for all provided times in
  * seconds with an initial composition defined by a map.
  *
  * \param chain: The chain to deplete with.
  * \param ccMap: A vector of initial concentration. Length of the vector
  * should match the length of the chain's nuclides attribute.
  * \param times: The length of the time step in seconds.
  * \param cutoff: Value under which a
  * concentration should be rounded to zero. default = 0.
  */
 std::vector<Eigen::VectorXd> run(const Chain &chain,
                                  const std::map<std::string, double> &ccMap,
                                  const std::vector<double> &times,
                                  double cutoff = 1.e-10);

 /**
  * \brief Solve the Bateman equation for all provided times in
  * seconds with an initial composition defined by a vector.
  *
  * \param chain: The chain to deplete with.
  * \param ccVector: A map of initial concentration in the ccMap["nuclide"] =
  * concentration format
  * \param times: The length of the time step in seconds.
  * \param cutoff: Value under which a concentration should be rounded to zero.
  * default = 0.
  */
 std::vector<Eigen::VectorXd> run(const Chain &chain, const Eigen::VectorXd &ccVector,
                                  std::vector<double> times,
                                  double cutoff = 1.e-10);

 Results results_;

private:
 /**
  * \brief Solve the Bateman equation on a timestep of size dt seconds with an
  * initial composition defined by a map.
  *
  * \param chain: The chain to deplete with.
  * \param ccVector: A vector of initial concentration. Length of the vector
  * should match the length of the chain's nuclides attribute.
  * \param dt: The length of the time step in seconds.
  * \param cutoff: Value under which a
  * concentration should be rounded to zero. default = 0.
  */
 std::vector<Eigen::VectorXd> run(const Chain &chain, const Eigen::VectorXd &ccVector,
                                  double dt, double cutoff = 1.e-10) const;

 /**
  * \brief theta_i tabulated values for CRAM-48
  *
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
  cdouble(+1.316284237125190e+1, +2.042951874827759e+1)
 };

 /**
  * \brief alpha_i tabulated values for CRAM-48
  *
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
  cdouble(+1.041366366475571e+2, -2.777743732451969e+2)
 };

 /**
  * \brief alpha_0 tabulated values for CRAM-48
  *
  */
 double alpha48_0 = 2.258038182743983e-47;
};

#endif
