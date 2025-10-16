from mpmath import mp, mpf, mpc, zeta, findroot, sin, log, fsum, linspace, diff
from sympy import nextprime
import time

mp.dps = 50  # Set precision to 50 decimal digits


def generate_primes(n):
    """Generate the first n prime numbers using sympy's nextprime."""
    primes = []
    current = mpf(1)
    for _ in range(n):
        current = nextprime(current)
        primes.append(current)
    return primes


def imaginary_sum(gamma, primes):
    """Compute the annihilation field A(gamma) = d sin(gamma * log p) for all primes."""
    return fsum([sin(gamma * log(p)) for p in primes])


def find_imaginary_zeros(gamma_range, primes, resolution=5000):
    """Find zero crossings of the annihilation field in the given range."""
    gamma_vals = linspace(mpf(gamma_range[0]), mpf(gamma_range[1]), resolution)
    im_vals = [imaginary_sum(g, primes) for g in gamma_vals]

    zero_crossings = []
    for i in range(len(gamma_vals) - 1):
        if im_vals[i] * im_vals[i + 1] < 0:
            try:
                root = findroot(
                    lambda g: imaginary_sum(mpf(g), primes),
                    (gamma_vals[i], gamma_vals[i + 1]),
                    solver='bisect',
                    verbose=False
                )
                zero_crossings.append(root)
            except Exception:
                continue
    return zero_crossings


def refine_zeta_zero(gamma_seed, tol=1e-30, maxsteps=10):
    """Refine gamma_seed using Newton's method on Re(zeta(1/2 + i*gamma))."""
    def f(g): return zeta(mpc(mpf('0.5'), g)).real
    def df(g): return diff(f, g)

    gamma = mpf(gamma_seed)
    for _ in range(maxsteps):
        try:
            # Added a try-except block to handle potential division by zero
            # if the derivative is flat, which can happen.
            f_gamma = f(gamma)
            df_gamma = df(gamma)
            if df_gamma == 0:
                return None
            delta = f_gamma / df_gamma
            gamma -= delta
            if abs(delta) < tol:
                return gamma
        except Exception:
            return None
    return None


def compute_zeta_zeros_in_range(
    gamma_range,
    num_primes,
    resolution,
    tol,
    outfile,
    start_idx
):
    """
    Main routine: find annihilation zeros, refine them, and classify true zeta zeros.
    MODIFIED: This function now appends to the file and returns the new zero count.
    """
    primes = generate_primes(num_primes) 
    seeds = find_imaginary_zeros(gamma_range, primes, resolution) 

    # This set will keep track of refined zeros found in this chunk to avoid duplicates
    # within the same run, using a string representation to handle precision issues.
    refined_in_session = set()
    idx = start_idx

    # MODIFICATION: Open the file in append mode ('a') to avoid overwriting previous results.
    with open(outfile, 'a') as data:
        for gamma in seeds:
            r = refine_zeta_zero(gamma)
            if r:
                # Check if this zero is a new discovery in this session
                r_str = str(r)
                if r_str in refined_in_session:
                    continue

                z = zeta(mpc(mpf('0.5'), r))
                # MODIFICATION: Relaxed duplicate check to only be within this session.
                # A proper duplicate check across sessions would require reading the file.
                if (abs(z) < tol and gamma_range[0] <= r <= gamma_range[1]):
                    refined_in_session.add(r_str)
                    msg = f"{idx}: Zeta Zero {r}"
                    print(msg)
                    data.write(msg + "\n")
                    idx += 1
    # Return the next index to start counting from
    return idx


if __name__ == "__main__":
    # --- Configuration ---
    START_GAMMA = 10
    CHUNK_SIZE = 100  # How large of a gamma range to search in each step
    NUM_PRIMES = 50
    RESOLUTION = 50000
    TOLERANCE = 1e-20
    OUTFILE = 'zeta_zeros.txt'

    current_gamma_start = mpf(START_GAMMA)
    zero_count = 1
    
    # Create or clear the output file at the beginning of the session
    with open(OUTFILE, 'w') as f:
        f.write(f"--- Zeta Zero Search Initialized at {time.ctime()} ---\n")

    print("--- Starting continuous search for Riemann zeta function zeros ---")
    print(f"Results will be saved to {OUTFILE}")
    print("Press Ctrl+C to stop the search gracefully.")

    # --- Main Loop ---
    try:
        while True:
            search_range = (current_gamma_start, current_gamma_start + CHUNK_SIZE)
            print(f"\n--- Searching in range: {search_range} ---")

            zero_count = compute_zeta_zeros_in_range(
                gamma_range=search_range,
                num_primes=NUM_PRIMES,
                resolution=RESOLUTION,
                tol=TOLERANCE,
                outfile=OUTFILE,
                start_idx=zero_count
            )

            # Prepare for the next chunk
            current_gamma_start += CHUNK_SIZE

    except KeyboardInterrupt:
        print(f"\n--- Search interrupted by user at {time.ctime()}. ---")
        print("Exiting gracefully.")
