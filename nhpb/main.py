from h.py import Lindbladian, Domain, evolve

def main():
    """
    calculate avergae photon number and second-order correlation
    nonHermitian Jaynes-Cummings Hamiltonian
    """
    # please input coupling constant and detuning parameters
    g = 0.005*67.5;  # [meV]
    det = 0;    # [meV]
    try:
        x = Domain();
        pars = Parameter(g, det);
        system = Lindbladian(pars);
        # system.update_pars(g, det);
        N, g2 = evolve(x, system);
        print(N, g2);
    except (TypeError, ValueError, ZeroDivisionError, AttributeError) as e:
        print(f"Error running calculation: {e}");
    except Exception as e:
        print(f"Unexpected error: {e}");

    return

if __name__ == "__main__":
    main();
