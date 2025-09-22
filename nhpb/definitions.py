# definitions
import math

from qutip import basis, tensor, destroy, qeye, mesolve, liouvillian, \
correlation_3op_2t

from . import par

class Par():
    """scaled parameters"""
    def __init__(self, detuning, coupling_const, kappa, gamma):
        """
        detuning -- detuning between cavity and spin mode [eV]
        coupling_const -- interaction strength between modes [< 0.9 eV]
        kappa -- T1 energy dissipative rate of cavity mode [eV]
        gamma -- T1 energy dissipative rate of spin mode [eV]
        """
        self.det, self.g, self.kappa, self.gamma = par.assign(
                detuning, coupling_const, kappa, gamma 
        )
        self.det_wc = par.freqs("cavity", self.det)
        self.det_ws = par.freqs("spin", self.det)

    def scaled(self, val):
        """scale paremeter value p wrt reference exchange frequency"""
        wr = max(self.g, self.kappa, self.gamma); # ref freq.
        return val/wr

    def times(self):
        N = 3*(2*math.pi); # number of cycles [rad]
        n = 90; # number of points per cycle
        tau = 1 / self.scaled(self.g)    # coherence (or observation) time length
        t_max = N * tau;
        w_max = max(
            self.scaled(self.gamma), self.scaled(self.kappa), 
            self.scaled(self.g), self.scaled(self.det)
        )
        T_min = 1 / w_max;
        dt = T_min / n;
        N_steps = int(t_max / dt);
        return [i*dt for i in range(N_steps + 1)]


class Operator():
    def __init__(self):
        # Precompute operators
        self._a = tensor(destroy(par.N), qeye(2));  # mode 1 (cavity)
        self._sm = tensor(qeye(par.N), destroy(2)); # mode 2 (qubit)

    def update_operator(self, mode1, mode2):
        self._a = mode1;
        self._sm = mode2;

    def Na(self): return self._a.dag() * self._a;  # occupation number operator

    def Nsm(self): return self._sm.dag() * self._sm;
    
    def X(self): 
        """coupling operator"""
        return self._a.dag()*self._sm + self._a*self._sm.dag()

    def D(self): return self._sm + self._sm.dag();  # displacement operator

    def JC_H(self, p):
        """
        Jaynes-Cummings Hamiltonian in RWA
        p -- class Par
        """
        H = [];
        for i in range(len(p.det_ws)):
            H.append( 
                p.scaled(p.det_wc[i]) * self.Na()
                    + p.scaled(p.det_ws[i]) * self.Nsm() 
                    + p.scaled(p.g) * self.X()
                    + p.scaled(par.v) * self.D()
            );
        return H

    def cops(self, p):
        """collapse operator with n_th = 0"""
        sqrt_kappa = p.scaled(p.kappa) ** 0.5;
        sqrt_gamma = p.scaled(p.gamma) ** 0.5;
        return [sqrt_kappa*self._a, sqrt_gamma*self._sm];


class Evolve():
    """Evolve density matrix in time"""
    def __init__(self, times, H, cops):
        """
        times -- list: floats
        H -- list: qobj
        cops -- list: qobj
        """
        self._psi0 = tensor(basis(par.N, par.nc), basis(2, par.nsm));
        self._t = times;    # time domain
        self._H = H
        self._cops = cops

    def occupation(self, Na, Nsm):
        """
        Na, Nsm -- qobj, occupation number operator
        """
        res = [];
        for w, _ in enumerate(self._H):
            # print("Hamiltonian at index", w, ":", H[w], type(H[w]))
            res.append(mesolve(self._H[w], self._psi0, self._t, self._cops,
                [Na, Nsm]));
        return res

    def correlation(self, op_1, op_2): 
        """
        two-time second-order correlation function
        <op_n(t)op_ndag(t+tau)op_n(t+tau)op_n(t)>
        Input
        op_n - qobj
        Output
        second order correlation at tau=0 -- list of list of type: qutip
        corr_mat : 1d-array, for op_n (first index) 
            at specified drive frequency (second index)
        """
        tau = [0]; # decoherence times
        a_ops = [op_1, op_2]
        b_ops = [op_1.dag()*op_1, op_2.dag()*op_2]

        corr_w = [];  # correlation versus drive freq
        corr_op = [];  # versus operator type
        for i, a_op in enumerate(a_ops):
            for w, _ in enumerate(self._H):
                corr_w.append(
                    correlation_3op_2t(
                        self._H[w], self._psi0, self._t, tau,
                        self._cops, a_op, b_ops[i], a_op
                    )[:,0]
                );
            corr_op.append(corr_w);
            corr_w = []
        return corr_op

############1234567
