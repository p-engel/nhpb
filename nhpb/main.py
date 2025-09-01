from .lindbladian import Lindbladian 

def main():
    """
    calculate avergae occupation number and second-order correlation
    nonHermitian Jaynes-Cummings Hamiltonian
    """
    try:
        L = Lindbladian();
        N = L.occupation()
        print(N[0].expect[0]);
    except (TypeError, ValueError, ZeroDivisionError, AttributeError) as e:
        print(f"Error running calculation: {e}");
        return None
    except Exception as e:
        print(f"Unexpected error: {e}");

    return N

if __name__ == "__main__":
    main();
