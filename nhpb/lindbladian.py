# delcaration in namespace of type Lindbladian

from .definitions import Operator, Evolve, Par

class Lindbladian(Operator):
    def __init__(self, detuning, coupling_const, kappa, gamma, drive, phase):
        """
        hamiltonian: (Qobj)
        cops: list of collapse operators (Qobj)
        """
        super().__init__()
        self.p = Par(detuning, coupling_const, kappa, gamma, drive, phase);  # parameter
        self._L = None;  # Cache for the superoperator

    def dynamics(self):
        """Evolve density matrix"""
        return Evolve(self.p.times(), self.JC_H(self.p), self.cops(self.p))
    
    def occupation(self):
        """average occupation number in time"""
        return self.dynamics().occupation(self.Na(), self.Nsm())

    def correlation(self):
        """compute second order correlation"""
        return self.dynamics().correlation(self._a, self._sm)
