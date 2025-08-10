# define functions

class Operator():
    par = Parameter();


    def __init__(self):
        """Hilbert space"""
        self.__basis1 = destroy(par.N);
        self.__basis2 = destroy(2);


    def update_Hilbert(basis1, basis2):
        self.__basis1 = basis1;
        self.__basis2 = basis2;


    def mode1(self):
        """cavity mode operator"""
        a = tensor(self.__basis1, qeye(2));
        
        return a


    def mode2(self):
        """qubit mode operator"""
        sm = tensor(qeye(par.N), self.__basis2);
        
        return sm


    def H(self, w):
        """
        JC Hamiltonian in RWA
        w - drive frequency
        """
        Na = self.mode1()*self.mode1().dag();   # number operator
        Nsm = self.mode2()*self.mode2().dag();
        X = self.mode1().dag()*self.mode2() + self.mode1()*self.mode2().dag();
        D = self.mode2() + self.mode2().dag();
        
        H = par.det_wa(w)*Na + par.det_ws(w)*Nsm + par.g*X + par.v*D;
        
        return H


    def cops(self):
        """collapse operator"""
        sqrt_kappa = par.kappa ** 0.5;
        sqrt_gamma = par.gamma ** 0.5;
        cops = [sqrt_kappa*self.mode1(), sqrt_gamma*self.mode2()];   # n_th = 0
        
        return cops
