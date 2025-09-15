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
v = 0
N = 12  # Size of Fock-space
nc = 0; nsm=1   # level in Hilbert space cavity and spin

def test(detuning: float, coupling_const: float, 
        kappa: float, gamma: float) -> float:
    try:
        res = w0 + detuning;
        res_cond = coupling_const < res;
        assert res_cond, "the model satisfies coupling constants smaller than \
            resonant frequency g < w < 1, see par file to change modal \
            frequencies."
        return detuning, coupling_const, kappa, gamma
    except AssertionError as e:
        print(f"Error: {e}")

def assign(detuning, coupling_const, loss_rate1, loss_rate2):
    if any(x is None for x in 
            [detuning, coupling_const, loss_rate1, loss_rate2]):
        return det, g, kappa, gamma
    else:
        test(detuning, coupling_const, loss_rate1, loss_rate2)

def freqs(mode, detuning):
    wlist = [w0, (w0+det)]
    match mode:
        case str("cavity"):
            return [ (w0 - w) for w in wlist]
        case str("spin"):
            return [ (w0+det - w) for w in wlist ]
        case _:
            return "Invalid mode"

## parameter based on measured Fano system in Nano Lett. by Masiello group
# ws = 0.9326
# det = 20e-3
# wc = ws - det
# gamma = 67.5e-3
# kappa = 9.126e-10
# g = 0.332e-3    # coupling constant
# v = 0.2 * 67.5e-3    # drive amplitude
# N = 12  # Size of Fock-space
# nc = 0; nsm=0   # level in Hilbert space cavity and spin
