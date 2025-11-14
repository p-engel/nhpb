# System parameters in eV

# constant
hbar = 1
radf_eV = 1.52e15   # [rad/s * 1/eV]

# test parameters
w0 = 1
det = 0  # detuning between spin and cavity modes
g = 0.05 # coupling constant
kappa = 0.005 / 2
gamma = 0.05 / 2
v = 0; phi = 0;
N = 12  # Size of Fock-space
nc = 0; nsm=0   # level in Hilbert space cavity and spin

def test(detuning: float, coupling_const: float, 
        kappa: float, gamma: float, drive: float, phase: float) -> float:
    try:
        res = w0 + detuning;
        res_cond = coupling_const < res;
        assert res_cond, "the model satisfies coupling constants smaller than \
            resonant frequency g < w < 1, see par file to change modal \
            frequencies."
        phase_cond = phase < 6.28
        assert phase_cond, "the phase is the range 0 to 2pi."
        return detuning, coupling_const, kappa, gamma, drive, phase
    except AssertionError as e:
        print(f"Error: {e}")

def assign(
    detuning, coupling_const, kappa, gamma, drive, phase
    ):
    if any(x is None for x in
            [detuning, coupling_const, kappa, gamma, drive, phase]):
        return det, g, kappa, gamma, v, phi
    else:
        return test(detuning, coupling_const, kappa, gamma, drive, phase)

def freqs(mode, detuning):
    wlist = [w0, (w0+det)]
    match mode:
        case str("cavity"):
            return [ (w0 - w) for w in wlist]
        case str("spin"):
            return [ (w0+det - w) for w in wlist ]
        case _:
            return "Invalid mode"

