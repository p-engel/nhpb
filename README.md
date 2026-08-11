# NHPB

Photon statistics of a weakly anharmonic mode coupled to a driven,
lossy cavity.

This repository contains the numerical and analytic approach
accompanying a manuscript in preparation. It computes steady-state
photon-number statistics of an anharmonic bosonic mode (e.g. a
transmon) coupled to a driven broadband bosonic (cavity) mode, in the
regime where the anharmonicity is small compared to the cavity
linewidth.

## Objectives

- Compute the steady state single- and two-excitation amplitudes of
  the coupled Kerr-cavity system under continuous-wave driving.

- Evaluate the second-order correlation function $g^{(2)}(0)$ of the
  light emitted through the cavity mode, and characterize the regime
  in which it becomes sub-Poissonian.

- Cross-check closed-form results derived under a few-excitation
  ansatz against a full Lindblad master-equation solution.

- Map the parameter regime over which the sub-Poissonian behaviour
  survives.

## Model

A bosonic mode `a` with Kerr nonlinearity,

```math
H_{NL} = \omega_a a^{\dagger}a + (\alpha/2) a^{\dagger 2}a^{2}
```

is coupled with strength $g$ to a cavity mode `c` driven by a CW laser
of amplitude $E$ and frequency $\omega$,

```math
\begin{aligned}
H_I &= g (a^{\dagger}c + c^{\dagger}a)
\\
H_c &= \omega_c c^{\dagger}c
  + (E e^{-i\omega t} + E^* e^{i\omega t})(c + c^{\dagger})
\end{aligned}
```

In the frame rotating at the laser frequency, with detunings
$\Delta_{a,c} = \omega_{a,c} - \omega$, the system Hamiltonian is

```math
H = \Delta_a a^{\dagger}a + (\alpha/2) a^{\dagger 2}a^{2}
  + \Delta_{c} c^{\dagger}c + g (a^{\dagger}c + c^{\dagger}a)
  + (E + E^*)(c + c^{\dagger})
```

The regime of interest is $\gamma \ll g < \kappa$, where $\gamma$
is the loss rate of the narrow discrete transmon mode coupled to a
broad continuum-like cavity mode, under weak drive $E \ll \kappa$,
with weak nonlinearity $\alpha \ll \kappa$.

## Formalism

Dissipation enters through a Lindblad master equation for the joint
density operator,

```math
\begin{aligned}
\dot{\rho} &= -i[H, \rho] + \gamma\mathcal{D}[a]\rho
  + \kappa\mathcal{D}[c]\rho
\\
\mathcal{D}[o] &= o \rho o^{\dagger} - \tfrac{1}{2}\{o^{\dagger} o, \rho\}
\end{aligned}
```

Steady states are obtained numerically with QuTiP. In parallel, an
effective non-Hermitian Hamiltonian (quantum-jump terms neglected)
together with a few-excitation ansatz truncated at two excitations,
yields closed-form expressions for the excitation amplitudes. The
numerics validate the closed-form expression.

## Observables

- $\langle c \rangle$ - steady-state single-photon amplitude in the
cavity mode

- $\langle c c \rangle$ - steady-state two-photon amplitude

- $N = \vert \langle c \rangle \vert^2$ - mean intracavity photon
  number

- $N / N_{g=0}$ - photon number normalised to the uncoupled
  ($g=0$) system

- $g^{(2)}(0)$ - equal-time second-order correlation, evaluated
  both from the ansatz amplitudes and from the normal-ordered Lindblad
  steady state

## Layout
 
```
nhpb/
  kerr_system.py       system construction: Hamiltonian, parameters, operators
  kerr_measurement.py  observables: amplitudes, photon number, g2
  lindbladian.py       master-equation assembly and steady-state solution
  par.py               parameter handling
  main.py              entry point
notebooks/
  analysis.ipynb       derivations, numerics, and manuscript figures
unit_test/             tests
```

## Installation
 
```bash
git clone https://github.com/p-engel/nhpb.git
cd nhpb
python3.12 -m venv .venv
source .venv/bin/activate
pip install -e .
```

Requires Python ≥ 3.10, QuTiP 5.x.

## Usage

```python
from nhpb.kerr_system import KerrSystem
from nhpb.kerr_measurement import KerrMeasurement
 
# all rates in units of the cavity linewidth κ_c
system = KerrSystem(
    delta_q=0.0,      # detuning of the nonlinear mode
    delta_cq=0.5,     # detuning between the two modes
    g=0.02,           # mode coupling
    kappa_c=1.0,      # cavity linewidth
    kappa_q=1e-4,     # linewidth of the nonlinear mode
    drive=1e-2,       # CW drive amplitude
    kerr=0.0,         # Kerr nonlinearity U
)
 
m = KerrMeasurement(system)
m.psi1()              # single-photon amplitude
m.psi2()              # two-photon amplitude
m.g2()                # g2(0), mean-field
m.g2_normal_ordered() # g2(0), normal-ordered Lindblad steady state
```

Parameter names in the code map to the notation above as follows:
 
| Code | Text | Meaning |
| --- | --- | --- |
| `kerr` | $\alpha$ | Kerr nonlinearity |
| `kappa_q` | $\gamma$ | linewidth of the anharmonic mode |
| `kappa_c` | $\kappa$ | cavity linewidth |
| `delta_q` | $\Delta_a$ | detuning of the anharmonic mode from the drive |
| `delta_cq` | $\delta_{ca}$ | detuning between the two modes |

## Status
 
Work in progress, accompanying a manuscript in preparation. Results
and derivations in `notebooks/analysis.ipynb` are subject to
revision. A citation and preprint link will be added on submission.

<!-- A mechanism describing non-Hermitian photon blockade

## Use Case
The system is described by a driven Jaynes-Cummings Hamiltonian. The TLS (spin
qubit) has a lifetime of ~10 femto-seconds, and the cavity, ~0.7 micro-seconds.

- Calculate the average occupation number in the qubit and cavity modes.
- Calculate the second-order autocorrelation through qubit and cavity modes.

-->