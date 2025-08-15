# definitions
import math

from qutip import tensor, destroy, qeye, mesolve, liouvillian 

import par

class Par():
    """scaled parameters"""
    def __init__(self, freq=None):
        """
        scale parameters by fastest dissipative rate value
        freq -- list: float, drive frequency
        """
        self.wr = max(par.kappa, par.gamma); # ref freq.
        self.gamma = par.gamma / self.wr
        self.kappa = par.kappa / self.wr
        self.g = par.g / self.wr
        self.v = par.v / self.wr
        self.det_wc = _scaled_detuning(par.wc)
        self.det_ws = _scaled_detuning(par.ws)

    def _scaled_detuning(freq0):
        return [ (freq0 - freq)/self.wr for freq in [par.ws, par.wc] ]


class Operator():
    def __init__(self):
        # Precompute operators
        self._a = tensor(destroy(par.N), qeye(2));  # mode 1 (cavity)
        self._sm = tensor(qeye(par.N), destroy(2)); # mode 2 (qubit)

    def update_operator(self, mode1, mode2):
        self._a = mode1;
        self._sm = mode2;

    def mode1(self): return self._a

    def mode2(self): return self._sm

    def occupation1(self): return self._a.dag() * self._a;

    def occupation2(self): return self._sm.dag() * self._sm;

    def JC_H(self, p):
        """
        Jaynes-Cummings Hamiltonian in RWA
        p -- class Par
        """
        coupling = self._a.dag()*self._sm + self._a*self._sm.dag();
        displace = self._sm + self._sm.dag();

        H = [];
        for i in range(p.ws):
        H.append( p.det_wc[i] * self.occupation1()
                    + p.det_ws[i] * self.occupation2() 
                    + p.g * coupling
                    + p.v * displace );
        return H

    def cops(self, p):
        """collapse operator with n_th = 0"""
        sqrt_kappa = p.kappa ** 0.5;
        sqrt_gamma = p.gamma ** 0.5;
        return [sqrt_kappa*self._a, sqrt_gamma*self._sm];


class Evolve():
    """Evolve density matrix in time"""
    def __init__(self, p)
        """
        p -- class Par, parameters
        hamiltonian -- class Operator, hamiltonian
        collaps_ops -- class Operator, collaps operator
        """
        self._psi0 = tensor(basis(par.N, par.na), basis(2, par.nsm));
        self._t = _times(p);

    def _times(p):
        N = 10; # number of oscillations/decay
        n = 20; # number of points per oscillations/decay
        tau_max = max(1/p.g, 1/p.gamma);  # longest life time
        t_max = N * tau;
        w_max = max(p.gamma, p.g, p.v, p.ws, p.wc)
        T_min = 2 * math.pi / w_max;
        dt = T_min / n;
        N_step = int(t_max / dt);
        return [i*dt for i in range(N_steps + 1)]

    def occupation(self, H, cops, Na, Nsm):
        """
        occupation number in modes
        Na, Nsm -- class, Operator, occupation number operator
        """
        res = [];
        for w in range(H):
            res.append(mesolve(H[w], psi0, self._t, cops, [Na, Nsm]));
        return res

    def correlation(self): return

############1234567
