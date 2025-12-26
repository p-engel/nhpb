import numpy as np

from nhpb.kerr_system import KerrSystem
from nhpb.kerr_measurement import KerrMeasurement

def test_amplitude():
    """
    unit test of single photon amplitude
    """
    n = 0.5;  # [kappa_c]
    p = 1e-3;  # ratio kappa_q : kappa_c
    lim = ( (n**2 + 1)/n * p )**0.5 * 1/2;   # [kappa_c]
    system1 = KerrSystem(delta_cq=n, g=0, kappa_c=1, kappa_q=p, drive=1e-2)
    system2 = KerrSystem(delta_cq=n, g=lim*5e1, kappa_c=1, kappa_q=p, drive=1e-2)
    M1 = KerrMeasurement(system1);
    M2 = KerrMeasurement(system2);
    a1 = abs(M1.psi1())**2;
    a2 = abs(M2.psi1())**2;
    a22 = abs(M2.psi2())**2
    try:
        expected_a1 = (1e-2)**2 / (n**2 + (1/2)**2);
        expected_a2 = 0;
        print(a1, expected_a1)
        print(a2, expected_a2)
        print(a22, a2**2)
        assert np.isclose(a1, expected_a1, rtol=1e-2), (
            f"the modulus square of the single photon amplitude should equal",
            "the occupation number"
            )
        assert np.isclose(a2, expected_a2, rtol=1e-2), (
            f"the modulus square of the steady state single photon amplitude",
            "should approach zero at the Fano antiresonance"
            )
        assert np.isclose(a22, a2**2, rtol=1e-2), (
            f"the modulus square of the steady state two-photon amplitude",
            "should equal the square of the single-photon amplitude",
            "for a linear system"
            )
    except AssertionError() as a:
        print(f"AssertionError : {a}")
    return

test_amplitude()
