import numpy as np
import h5py


class SolverOrdered:
    # def __init__(self, matrix, isotopes: list, cc: dict):
    def __init__(self, chain):
        self.chain = chain
        self.isotopes, self.matrix = chain.build_decay_matrix("topo")
        self.F = None
        self.Ns = None
        self.cc = {}
        self.source = {}
        self.func = {}

    def set_initial_cc(self, cc: dict):
        self.cc = cc

    def set_source(self, source: dict):
        self.source = source

    def compute_coeffs(self, verbose=False):
        # Checking that the matrix is lower triangular
        X, Y = self.matrix.nonzero()
        assert np.all(X >= Y)

        isotopes = [(i, val) for i, val in enumerate(self.isotopes)]
        Ns = {key: {} for key in self.isotopes}
        F = {key: {} for key in self.isotopes}

        for i, i_isot in isotopes:
            stable = self.chain.stable(i_isot)
            _, Y = self.matrix[[i]].nonzero()
            Y = Y[:-1]
            parents = [isotopes[y] for y in Y]

            # Computing Fik
            for k, k_isot in isotopes[:i][::-1]:  # 0 <= k < i
                factor = - 1 / (self.matrix[i, i] - self.matrix[k, k])
                Fikm = 0
                Fikp = 0
                for j, j_isot in isotopes[k:i]:
                    if (Fiktemp := self.matrix[i, j] * F[j_isot][k_isot]) >= 0:
                        Fikp += Fiktemp
                    else:
                        Fikm += Fiktemp
                Fik = factor * (Fikm + Fikp)
                F[i_isot][k_isot] = Fik
                # F[i_isot][i_isot] -= Fik

            # Computing Ns(i)
            if stable:
                Ns[i_isot] = self.cc.get(i_isot, 0) - sum([F[i_isot][k_isot] for _, k_isot in isotopes[:i]])
            else:
                Ns[i_isot] = - sum([self.matrix[i, k] * Ns[k_isot] for k, k_isot in parents]) / self.matrix[i, i]
                Ns[i_isot] -= self.source.get(i_isot, 0.) / self.matrix[i, i]

            # Computing Fii
            if stable:
                F[i_isot][i_isot] = sum([self.matrix[i, j] * Ns[j_isot] for j, j_isot in parents])
                F[i_isot][i_isot] += self.source.get(i_isot, 0.)
            else:
                F[i_isot][i_isot] = self.cc.get(i_isot, 0)
                F[i_isot][i_isot] -= Ns[i_isot]
                F[i_isot][i_isot] -= sum([F[i_isot][k_isot] for k, k_isot in isotopes[:i][::-1]])

            # Building the resulting function
            if stable:
                self.func[i_isot] = (lambda t, isot_=i_isot: Ns[isot_]
                                     + sum([fik_ * np.exp(-self.chain.dconst(k_) * t) for k_, fik_ in self.F[isot_].items() if k_ != isot_])
                                     + self.F[isot_][isot_] * t
                                     )
            else:
                self.func[i_isot] = (lambda t, isot_=i_isot: Ns[isot_]
                                     + sum([fik_ * np.exp(-self.chain.dconst(k_) * t) for k_, fik_ in self.F[isot_].items()]))

        self.F = F
        self.Ns = Ns

    def write_coeffs(self, path):
        with h5py.File(path, 'w') as f:
            for i_isot in self.F:
                for j_isot in self.F[i_isot]:
                    if self.F[i_isot][j_isot] != 0:
                        f[f"Fij/{i_isot}/{j_isot}"] = np.array([self.F[i_isot][j_isot]])

        #     for kid_back, k in enumerate(Y[::-1]):
        #         k_id = len(Y) - kid_back - 1
        #         parent = self.isotopes[k]
        #         if verbose:
        #             print(f"\t{parent=}")
        #         if i != k:
        #             Cik = self.matrix[i, k]
        #             C[i_isot][parent] = Cik
        #             Ns[i_isot] += Cik * Ns[parent]
        #             denom = C[i_isot][i_isot] - C[parent][parent]
        #             Fikm = 0
        #             Fikp = 0
        #             if verbose:
        #                 print("\tC(%s, %s) = %.8e" % (i_isot, parent, Cik))
        #                 print("\tDenom = %.4e - %.8e = %.e" % (C[i_isot][i_isot], C[parent][parent],denom))
        #             for j_id, j in enumerate(Y[k_id:-1]):
        #                 jisot = self.isotopes[j]
        #                 try:
        #                     Fik_ = C[i_isot][jisot] * F[jisot][parent]
        #                     if Fik_ < 0:
        #                         Fikm += Fik_
        #                     else:
        #                         Fikp += Fik_
        #                     if verbose:
        #                         print("\t\tF[%s][%s] = %8e" % (jisot, parent, F[jisot][parent]))
        #                 except KeyError:
        #                     pass
        #             Fik = (Fikm + Fikp) / denom
        #             if verbose:
        #                 print("\tF[%s][%s] = %.8e" % (i_isot, parent, Fik))
        #             F[i_isot][parent] = Fik

        #     Ns[i_isot] /= C[i_isot][i_isot]

        #     F[i_isot][i_isot] = self.cc.get(i_isot, 0) - Ns[i_isot] 
        #     for j_id, j in enumerate(Y[:-1]):
        #         jisot = self.isotopes[j]
        #         F[i_isot][i_isot] -= F[i_isot][jisot]
        #     if verbose:
        #         print("\tF[%s][%s] = %.8e" % (i_isot, i_isot, F[i_isot][i_isot]))
        # self.F = F
        # self.C = C
        # self.Ns = Ns
