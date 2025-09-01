# delcaration in namespace of type Lindbladian

from definitions import Operator, Evolve, Par

class Lindbladian(Operator):
    def __init__(self):
        """
        hamiltonian: (Qobj)
        cops: list of collapse operators (Qobj)
        """
        super().__init__()
        self.p = Par();  # parameter
        self._L = None;  # Cache for the superoperator
    
    def build(self):
        """Construct the Lindbladian superoperator"""
        self._L = liouvillian(self.JC_H(self.p), self.cops(self.p));
        return self._L
    
    def apply(self, rho):
        """Apply L[rho]"""
        if self._L is None:
            self.build();
        return self._L @ rho.full().ravel()  # returns vectorized form

    def dynamics(self):
        """Evolve density matrix"""
        return Evolve(self.p.times())
    
    def occupation(self):
        """average occupation number in time"""
        return self.dynamics().occupation(self.JC_H(self.p), self.cops(self.p),
                                            self.occupation1(),
                                            self.occupation2())

    def correlation(self):
        """compute second order correlation"""
        return self.dynamics().correlation()
