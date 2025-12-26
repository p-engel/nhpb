import numpy as np
from qutip import basis, tensor, destroy, qeye

# ============================================================
# Class — KerrSystem
# ============================================================

class KerrSystem:
    """
    Represents a driven cavity coupled to a Kerr nonlinear cavity.
    Defines operators, Hamiltonian components, and collapse operators.
    """

    def __init__(
        self,
        N_c: int = 10,
        N_q: int = 10,
        delta_cq: float = 0,
        delta_q: float = 0,
        g: float = 0.0,
        kappa_c: float = 0.0,
        kappa_q: float = 0.0,
        drive: float = 0.0,
        kerr: float = 0.0,
    ):
        """
        Parameters [units of rate]
        ----------
        N_c         : cavity Hilbert dimension
        N_q         : Kerr cavity Hilbert dimension
        delta_cq    : linear--nonlinear cavity detuning frequency
        delta_q     : nonlinear cavity--drive detuning frequency
        g           : coupling constant rate
        kappa_c     : cavity decay rate
        kappa_q     : qubit decay rate
        drive       : driving rate on cavity
        kerr        : Kerr nonlinearity strength
        """

        # store parameters
        self.N_c = N_c
        self.N_q = N_q
        self.delta_cq = delta_cq
        self.delta_q = delta_q
        self.g = g
        self.kappa_c = kappa_c
        self.kappa_q = kappa_q
        self.drive = drive
        self.kerr = kerr

        # ---------------------
        # Operators
        # ---------------------
        self.a_c = tensor(destroy(N_c), qeye(N_q))
        self.a_q = tensor(qeye(N_c), destroy(N_q))

        self.n_c = self.a_c.dag() * self.a_c
        self.n_q = self.a_q.dag() * self.a_q

    # ---------------------------
    # Hamiltonian pieces
    # ---------------------------

    def H_c(self): return self.delta_cq * self.n_c

    def H_qc(self):
        return self.delta_q * (self.n_q + self.n_c)

    def H_interaction(self):
        return self.g * (self.a_c * self.a_q.dag() +
                         self.a_c.dag() * self.a_q)

    def H_drive(self):
        """Drive term: ε (a e^{-iφ} + a† e^{iφ})"""
        return self.drive * (self.a_c + self.a_c.dag())

    def H_kerr(self):
        """Kerr term: χ a† a† a a"""
        return self.kerr * (self.a_q.dag() * self.a_q.dag() * 
                            self.a_q * self.a_q)

    def H(self):
        """Full system Hamiltonian"""
        return (
            self.H_c()
            + self.H_qc()
            + self.H_interaction()
            + self.H_drive()
            + self.H_kerr()
        )

    # -----------------------------
    # Collapse operators
    # -----------------------------

    def collapse_ops(self):
        cops = []
        if self.kappa_c > 0:
            cops.append(np.sqrt(self.kappa_c) * self.a_c)
        if self.kappa_q > 0:
            cops.append(np.sqrt(self.kappa_q) * self.a_q)
        return cops









































