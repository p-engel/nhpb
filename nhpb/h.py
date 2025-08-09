# Definitions
from def.py import Parameter, Domain

class Lindbladian():
	"""
	Struct for Lindbladian
    """
    def __init__(self, par):
        """
        default constructor, with private data as self
        """
        self.p = par;  [array, double]

    def psi0(self)
        """
        intial state
        """
        psi0 = tensor(basis(self.p.N, self.p.na), \
                    basis(2, self.par.nsm));

        return psi0

    def operator(self, w)
        """
        JC Hamiltonian in RWA
        w - drive frequency
        """
        a = tensor(destroy(self.p.N), qeye(2));
        sm = tensor(qeye(self.p.N), destroy(2));
        Na = a.dag()*a;     # number operator for cavity mode
        Ns = sm.dag()*sm;
        X = a.dag()*sm + a*sm.dag();    # coupling operator
        D = sm + sm.dag();  # displacement operator
        det_wa = self.p.wa - w;     # detuned freq. wrt drive
        det_ws = self.p.ws - w;
        g = self.p.g;
        v = self.p.E;
        sqrt_kappa = self.p.kappa ** 0.5;
        sqrt_gamma = self.p.gamma ** 0.5;
        
        H = det_wa*Na + det_ws*Ns + g*X + v*D;
        cops = [sqrt_kappa*a, sqrt_gamma*sm];   # n_th = 0

        return H, cops, a, sm

    def update_par(self, g, det)
        """
        modify parameters
        input - double
            g - coupling constant betweem modes
            det - detuning between modes
        """
        self.p = Parameter(coupling_const==g, detuning=det);

        return

####    output = mesolve(H, psi0, tlist, c_ops, [Na, Nsm])
