# Definitions
from def.py import parameters

class Lindbladian():
	"""
	Struct for Lindbladian
    """
    def __init__(self):
        """
        default constructor, with private data as self
        """
        self.par = parameters(detuning, g);  [array, double]
        self.w = domain().w; [array, double]

    def psi0(self)
        """
        intial state
        """
        psi0 = tensor(basis(self.par.N, self.par.na), \
                    basis(2, self.par.nsm))

        return psi0

    def operators(self)
        """
        JC Hamiltonian in RWA
        """
        a = tensor(destroy(self.par.N), qeye(2))
        sm = tensor(qeye(self.par.N), destroy(2))
        H = (self.par.wc - self.w)*a.dag()*a + \
                (self.par.wa - self.w)*sm.dag()*sm + \
                self.par.g*(a.dag()*sm + a*sm.dag()) + \
                self.par.E*(sm + sm.dag())
        c_ops = [np.sqrt(kappa)*a, np.sqrt(gamma)*sm]   # n_th = 0

        return H, cops, a, sm

    def update_par(self, det, g)
        """
        modify parameters
        input - double
            det - detuning between modes
            g - coupling constant betweem modes
        """
        self.par = parameters(detuning=det, coupling_const==g);

        return

####    output = mesolve(H, psi0, tlist, c_ops, [Na, Nsm])
