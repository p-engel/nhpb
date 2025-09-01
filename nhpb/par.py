# System parameters in eV
# based on measured Fano system in Nano Lett. by Masiello group

# constant
hbar = 1
radf_eV = 1.52e15   # [rad/s * 1/eV]

# parameter
ws = 0.9326
det = 20e-3
wc = ws - det
gamma = 67.5e-3
kappa = 9.126e-10
g = 0.332e-3    # coupling constant
v = 0.2 * 67.5e-3    # drive amplitude
N = 12  # Size of Fock-space
nc = 0; nsm=0   # level in Hilbert space cavity and spin

def detuning(freq0):
    return [ (freq0 - freq) for freq in [ws, wc] ]
