from lindbladian import Lindbladian 

def main():
    """
    calculate avergae occupation number and second-order correlation
    nonHermitian Jaynes-Cummings Hamiltonian
    """
    try:
        L = Lindbladian();
        N = L.occupation()
        print(N.expect[0], N.expect[1]);
    except (TypeError, ValueError, ZeroDivisionError, AttributeError) as e:
        print(f"Error running calculation: {e}");
    except Exception as e:
        print(f"Unexpected error: {e}");

    return

if __name__ == "__main__":
    main();
