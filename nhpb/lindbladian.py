# delcaration in namespace of type Lindbladian

from qutip import 

from def.py import Operator, Evolve, Domain

class Lindbladian(Operator):
    def __init__(self, hamiltonian=None, collapse_ops=None):
        """
        hamiltonian: (Qobj)
        state: initial state (Qobj)
        cops: list of collapse operators (Qobj)
        dim: optional, system dimension for validation
        """
        self._H = hamiltonian if hamiltonian is not None else self.JC_H(w);
        self._cops = collapse_ops if collapse_ops is not None else self.cops();
        self._dim = hamiltonian.shape[0] or Operator.JC_H()[0];
        self._L = None;  # Cache for the superoperator
    
    def build(self):
        """Construct the Lindbladian superoperator"""
        self._L = liouvillian(self._H, self._cops);
        return self._L
    
    def add_collapse(self, C):
        """Add a collapse operator"""
        self.cops.append(C);
        self._L = None;  # invalidate cached L
    
    def set_hamiltonian(self, hamiltonian):
        """Update the Hamiltonian"""
        self._H = hamiltonian;
        self._L = None;
    
    def apply(self, rho):
        """Apply L[rho]"""
        if self._L is None:
            self.build();
        return self._L @ rho.full().ravel()  # returns vectorized form
    
    def dynamics(self):
        """Evolve density matrix in time"""
        return Evolve(Domain())
