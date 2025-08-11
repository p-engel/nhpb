# definitions

from qutip import tensor, destroy, qeye, mesolve, liouvillian 

import par

class Operator():
    def __init__(self):
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

    def occupation1(self):
        """number operator in mode 1"""
        return self._a.dag() * self._a;

    def occupation2(self):
        return self._sm.dag() * self._sm;

    def JC_H(self, w):
        """Jaynes-Cummings Hamiltonian in RWA"""
        coupling = self._a.dag()*self._sm + self._a*self._sm.dag();
        displace = self._sm + self._sm.dag();
        return (  par.det_wa(w) * self.occupation1()
                + par.det_ws(w) * self.occupation2() 
                + par.g * coupling
                + par.v * displace  );

    def cops(self):
        """collapse operator with n_th = 0"""
        sqrt_kappa = par.kappa ** 0.5;
        sqrt_gamma = par.gamma ** 0.5;
        return [sqrt_kappa*self._a, sqrt_gamma*self._sm];


class Evolve():
    """Evolve density matrix in time"""
    def __init__(self, domain)
        """
        occupation - number operator in mode 1
        domain - obj-type with frequency and time list
        """
        self._psi0 = tensor(basis(par.N, par.na), basis(2, par.nsm));
        self._dom = domain;

    def occupation(self, Na, Nsm):
        """occupation number in modes"""
        res = [];
        for w in self.dom.w():
            H = self.op.JC_H(w);
            res.append(mesolve(H, psi0, self.dom.t(), self._cops, [Na, Nsm]));
        return res

    def correlation(self):
        return

############1234567
