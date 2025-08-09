from h.py import Lindbladian

def main():
    """
    calculate avergae photon number and second-order correlation
    nonHermitian Jaynes-Cummings Hamiltonian
    """
    # please input coupling constant and detuning parameters
    g = 0.005*67.5;  # [meV]
    det = 0;    # [meV]
    try:
        system = Lindbladian();
        system.update_pars(g, det);
        N = system.occupation();    # average occupation number
        g2 = system.correlation();    # second order correlation
        print(N, g2);
    except (TypeError, ValueError, ZeroDivisionError, AttributeError) as e:
        print(f"Error running calculation: {e}");
    except Exception as e:
        print(f"Unexpected error: {e}");

    return

if __name__ == "__main__":
    main();
