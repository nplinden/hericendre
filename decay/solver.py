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
