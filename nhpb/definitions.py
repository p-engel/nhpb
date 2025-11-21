# definitions
import math, cmath

from qutip import basis, tensor, destroy, qeye, mesolve, liouvillian, \
correlation_3op_2t, steadystate, expect

from . import par

# constants
cavity = 1
qubit = 2


class Par():
    """scaled parameters"""
    def __init__(self, detuning, coupling_const, kappa, gamma, drive, phase):
        """
        detuning -- detuning between cavity and spin mode [eV]
        coupling_const -- interaction strength between modes [< 0.9 eV]
        kappa -- T1 energy dissipative rate of cavity mode [eV]
        gamma -- T1 energy dissipative rate of spin mode [eV]
        """
        self.det, self.g, self.kappa, self.gamma, self.v, self.phi = par.assign(
                detuning, coupling_const, kappa, gamma, drive, phase
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
        # tau = 1 / self.scaled(self.g)    # coherence (or observation) time length
        tau = 1 / self.scaled(self.kappa)
        t_max = N * tau;
        w_max = max(
            self.scaled(self.gamma), self.scaled(self.kappa), 
            self.scaled(self.g), self.scaled(self.det), 
            self.scaled(self.v)
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

    def Dsm(self): return self._sm + self._sm.dag();  # displacement operator

    def Da(self, phi): 
        return self._a*cmath.exp(-1j*phi) + cmath.exp(1j*phi)*self._a.dag();

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
                    + p.scaled(p.v) * self.Da(p.phi)
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


class Kerr_sys():
    """The Kerr-Hamiltonian system"""
    # constants
    cavity = 1
    qubit = 2
    interaction = 3
    displace = 4
    kerr_nonlinear = 5

    def __init__(self, delta_q, delta_cq, u, g, v, kappa_q, kappa_c):
        self.wq = np.array(delta_q)
        self.det = delta_cq
        self.u = u
        self.g = g
        self.v = v
        self.kappa_q = kappa_q
        self.kappa_c = kappa_c
        self.ac = tensor(destroy(par.N), qeye(par.N));
        self.aq = tensor(qeye(par.N), destroy(par.N))
    
    def N(self, mode): 
        match mode
            case cavity: return self.ac.dag() * self.ac
            case qubit: return self.aq.dag() * self.aq
            case _: "Undefined mode"

    def X(self): 
        """coupling operator"""
        return self.ac.dag()*self.aq + self.aq.dag()*self.ac

    def D(self):
        """displacement operator"""
        return self.aq + self.aq.dag();

    def N2(self):
        """Kerr nonlinear operator"""
        return self.aq.dag() * self.aq.dag() * self.aq * self.aq

    def H(self, mode):
        match mode
            case cavity: return self.det * self.N(cavity)
            case qubit: return self.wq * self.N(qubit)
            case interaction: return self.g * self.X()
            case displace: return self.v * self.D()
            case kerr_nonlinear: return self.u/2 * self.N2()
            case _: "undefined Hamiltonian term"

    def Hs(self):
        """system Hamiltonian"""
        sys = self.H(cavity) + self.H(qubit) + self.H(interaction) \
            + self.H(kerr_nonlinear) + self.H(displace)
        return sys

    def cops(self):
        """collapse operator with n_th = 0"""
        sqrt_kappa_c = self.kappa_c ** 0.5;
        sqrt_kappa_q = self.Kappa_q ** 0.5;
        return [sqrt_kappa_c*self.ac, sqrt_kappa_q*self.aq];


class Evolve_ss():
    def __init__(self, H, cops):
        """steady state solution of the system's Lindbladian"""
        self.H = H
        self.cops = cops

    def rho_ss(self): return steadystate(self.H, self.cops)

    def expect_ss(self, op):
        """expectation value of operator <op>"""
        rho_f = self.rho_ss()
        return expect(op, rho_f)
        

############1234567
