# declaration for system of two-coupled linear and nonlinear bossons. 
# the linear mode is driven and the nonlinear mode has a Kerr-nonlinearity.
# the output mode is the linear mode's operations in steady state
from .definitions import Kerr_sys, Evolve_ss


class Kerr_sys():
    """The Kerr-Hamiltonian"""
    def __init__(self, wq, det, u, g, v, kappa_q, kappa_c):
        """
        hbar = 1
        wq - float: resonance energy of nonlinear mode
        det - float: detuning energy of linear and nonlinear mode
        u - float: energy proportional to Kerr nonlinearity
        g - float: interaction energy between the modes
        v - float: intercation energy between linear mode and CW laser
        kappa - float: dissipation energy of nonlinear and linear mode q, c
        """
        self.kerr_H = Kerr_sys(wq, det, u, g, v, kappa_q, kappa_c)

    def one_photon(self):
        """single photon amplitude in output mode"""
        return self.kerr_H.one_photon()  # [np 1d array]

    def two_photon(self):
        """two photon amplitude in output mode"""
        return self.kerr_H.two_photon()  # [np 1d array]

    def photon_number(self):
        """average photon number in output mode"""
        return self.kerr_H.N()  # [np 1d array]

    def g2():
        """normalized second order correlatoion"""
        mod_twophoton = np.absolute(self.two_photon)
        mod_onephoton = np.absolute(self.one_photon)
        g2 = (mod_twophoton ** 2) / (mod_onephoton ** 4)
        return g2  # [np 1d array]
