# define functions

class Operator():
    def __init__(self, params):
        """
        params: object holding system parameters
        """
        self.par = params;
        # Precompute operators
        self._a = tensor(destroy(self.par.N), qeye(2));  # mode 1 (cavity)
        self._sm = tensor(qeye(self.par.N), destroy(2)); # mode 2 (qubit)

    def update_operator(self, mode1, mode2):
        self._a = mode1;
        self._sm = mode2;

    def mode1(self):
        return self._a

    def mode2(self):
        return self._sm

    def JC_H(self, w):
        """Jaynes-Cummings Hamiltonian in RWA"""
        Na = self._a * self._a.dag();   # number operator
        Nsm = self._sm * self._sm.dag();
        X = self._a.dag()*self._sm + self._a*self._sm.dag();
        D = self._sm + self._sm.dag();
        return (  self.par.det_wa(w) * Na
                + self.par.det_ws(w) * Nsm 
                + self.par.g * X 
                + self.par.v * D  );

    def cops(self):
        """collapse operator with n_th = 0"""
        sqrt_kappa = self.par.kappa ** 0.5;
        sqrt_gamma = self.par.gamma ** 0.5;
        return [sqrt_kappa*self._a, sqrt_gamma*self._sm];

def psi0(par)
    """intial state"""
    return tensor(basis(par.N, par.na), basis(2, par.nsm));
