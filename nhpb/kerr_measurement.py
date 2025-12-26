from qutip import steadystate, expect
from .kerr_system import KerrSystem
# =================================================================
# KerrMeasurement — steady state observables in linear cavity mode
# =================================================================
class KerrMeasurement:
    """
    Given a KerrSystem, computes steady state, expectation values,
    and correlations g^(2).
    """
    def __init__(self, system: KerrSystem):
        self.sys = system
        self._rho_ss = None   # cached after first computation

    # --------------------
    # Core
    # --------------------
    def steady_state(self):
        """Compute and store steady-state density matrix."""
        if self._rho_ss is None:
            H = self.sys.H()
            cops = self.sys.collapse_ops()
            self._rho_ss = steadystate(H, cops)
        return self._rho_ss

    def expect_ss(self, op):
        """Expectation value in steady state."""
        return expect(op, self.steady_state())

    def psi1(self):
        """linear cavity single-photon amplitude <a>"""
        return self.expect_ss(self.sys.a_c)

    def psi2(self):
        """linear cavity two-photon amplitude <a a>"""
        return self.expect_ss(self.sys.a_c * self.sys.a_c)

    # ---------------------------------------------------------
    # Observable: occupation number, g^(2) correlation
    # ---------------------------------------------------------
    def occupation(self):
        """<a† a>"""
        return self.expect_ss(self.sys.n_cav)

    def g2(self):
        """g^(2) = |<a a>|^2 / |<a>|^4"""
        a1 = self.psi1()
        a2 = self.psi2()
        return abs(a2)**2 / abs(a1)**4

    def g2_normal_ordered(self):
        """<a† a† a a> / <a† a>^2"""
        a = self.sys.a_cav
        num = self.expect_ss(a.dag() * a.dag() * a * a)
        den = self.occupation() ** 2
        return num / den

