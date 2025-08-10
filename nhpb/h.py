# Definitions
from def.py import Parameter, Domain

class Lindbladian():
	"""
	Struct for Lindbladian
    """
    def __init__(self, par):
        self.p = par;  [array, double]

    def psi0(self)
        """
        intial state
        """
        psi0 = tensor(basis(self.p.N, self.p.na), \
                    basis(2, self.par.nsm));

        return psi0

####    output = mesolve(H, psi0, tlist, c_ops, [Na, Nsm])
