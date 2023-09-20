import numpy as np
import yaml
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from decay.utils import ROOT_PATH


class CramSolver:
    def __init__(self, chain, order=48):
        self.isotopes, self.matrix = chain.build_decay_matrix(order="zae")

        with open(f"{ROOT_PATH}/cram.yml") as f:
            coeffs = yaml.safe_load(f)[order]

        self.theta = np.array([r + i * 1j for r, i in zip(coeffs["theta_real"], coeffs["theta_imag"])], 
                              dtype=np.complex128)
        self.alpha = np.array([r + i * 1j for r, i in zip(coeffs["alpha_real"], coeffs["alpha_imag"])], 
                              dtype=np.complex128)
        self.alpha_0 = coeffs["alpha_0"]

    def single_step(self, n0, dt):
        matrix = sp.csr_matrix(self.matrix * dt, dtype=np.float64)
        y = n0.copy()
        identity = sp.eye(matrix.shape[0])
        for alpha, theta in zip(self.alpha, self.theta):
            y += 2*np.real(alpha*sla.spsolve(matrix - theta*identity, y))
        return y * self.alpha_0

    def __call__(self, n0: np.ndarray, times: np.ndarray):
        n = np.zeros(shape=(len(times), len(n0)))
        n[0] = n0
        dts = times[1:] - times[:-1]
        for idt, dt in enumerate(dts):
            n[idt+1] = self.single_step(n[idt], dt)
        return n

