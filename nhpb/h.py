# Definitions
import numpy as np
from qutip import Qobj, liouvillian

from def.py import Parameter, Domain

class Lindbladian:
    def __init__(self, H, psi0, collapse_ops=None, dim=None):
        """
        H: Hamiltonian (Qobj)
        psi0: initial state
        collapse_ops: list of collapse operators (Qobj)
        dim: optional, system dimension for validation
        """
        self.H = H
        self.psi0 = psi0
        self.collapse_ops = collapse_ops if collapse_ops is not None else []
        self.dim = dim or H.shape[0]
        self._L = None  # Cache for the superoperator
    
    def build(self):
        """Construct the Lindbladian superoperator."""
        self._L = liouvillian(self.H, self.collapse_ops)
        return self._L
    
    def add_collapse(self, C):
        """Add a collapse operator."""
        self.collapse_ops.append(C)
        self._L = None  # invalidate cached L
    
    def set_hamiltonian(self, H):
        """Update the Hamiltonian."""
        self.H = H
        self._L = None
    
    def apply(self, rho):
        """Apply L[rho]."""
        if self._L is None:
            self.build()
        return self._L @ rho.full().ravel()  # returns vectorized form
    
    def evolve(self, psi0, tlist, Na, Nsm):
        """Evolve density matrix in time (optional)."""
        from qutip import mesolve
        return mesolve(self.H, self.psi0, tlist, self.collapse_ops, [Na, Nsm])
