from .lindbladian import Lindbladian 

def main(detuning=None, coupling_const=None, kappa=None, gamma=None):
    """
    calculate avergae occupation number and second-order correlation
    nonHermitian Jaynes-Cummings Hamiltonian
    """
    try:
        L = Lindbladian(detuning, coupling_const, kappa, gamma);
        N = L.occupation()
        G2 = L.correlation()
        print(N[0].expect[0]);
        print(G2[0][0], type(G2[0][1]), G2[0][1].shape, 
            G2[0][0][-1], G2[0][1][-1]);
    except (TypeError, ValueError, ZeroDivisionError, AttributeError) as e:
        print(f"Error running calculation: {e}");
        return None
    except Exception as e:
        print(f"Unexpected error: {e}");

    return N

if __name__ == "__main__":
    main();
