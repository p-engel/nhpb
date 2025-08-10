# definitions
import par

class Operator():
    def __init__(self, par):
        # Precompute operators
        self._a = tensor(destroy(par.N), qeye(2));  # mode 1 (cavity)
        self._sm = tensor(qeye(par.N), destroy(2)); # mode 2 (qubit)

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
        return (  par.det_wa(w) * Na
                + par.det_ws(w) * Nsm 
                + par.g * X 
                + par.v * D  );

    def cops(self):
        """collapse operator with n_th = 0"""
        sqrt_kappa = par.kappa ** 0.5;
        sqrt_gamma = par.gamma ** 0.5;
        return [sqrt_kappa*self._a, sqrt_gamma*self._sm];


class Evolve():
    """Evolve density matrix in time"""
    def __init__(self)
    self._psi0 = tensor(basis(par.N, par.na), basis(2, par.nsm));
    self._op = Operator();
    self._dom = Domain();

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
