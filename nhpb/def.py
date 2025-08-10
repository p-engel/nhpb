# definitions

class Operator():
    def __init__(self, pars):
        # Precompute operators
        self._a = tensor(destroy(pars.N), qeye(2));  # mode 1 (cavity)
        self._sm = tensor(qeye(pars.N), destroy(2)); # mode 2 (qubit)

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
        return (  pars.det_wa(w) * Na
                + pars.det_ws(w) * Nsm 
                + pars.g * X 
                + pars.v * D  );

    def cops(self):
        """collapse operator with n_th = 0"""
        sqrt_kappa = pars.kappa ** 0.5;
        sqrt_gamma = pars.gamma ** 0.5;
        return [sqrt_kappa*self._a, sqrt_gamma*self._sm];

def state(pars):
    """intial state"""
    return tensor(basis(pars.N, pars.na), basis(2, pars.nsm));

class Evolve():
    """Evolve density matrix in time"""
    def __init__(self)
    self._psi0 = state(pars);
    self._op = Operator(pars);
    self._dom = Domain(pars);

    def occupation():
        """occupation number in modes"""
        Na = self.op.mode1() * self.op.mode1().dag();
        Nsm = self.op.mode2() * self.op.mode2().dag();            
        res = [];
        for w in self.dom.w():
            H = op.JC_H(w);
            res.append(mesolve(H, psi0, self.dom.t(), self._cops, [Na, Nsm]));
        return res

    def correlation():
        return

############1234567
