from mpmath import mp, mpf, mpc, zeta, findroot, sin, log, fsum, linspace, diff
from sympy import nextprime
from sys import argv

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
    """Refine gamma_seed using Newton's method on Re(?(1/2 + i?))."""
    def f(g): return zeta(mpc(mpf('0.5'), g)).real
    def df(g): return diff(f, g)

    gamma = mpf(gamma_seed)
    for _ in range(maxsteps):
        delta = f(gamma) / df(gamma)
        gamma -= delta
        if abs(delta) < tol:
            return gamma
    return None


def compute_zeta_zeros_with_classification(
    gamma_range=(10, 20),
    num_primes=100,
    resolution=5000,
    tol=1e-10,
    outfile='zeta_zeros.txt'
):
    """Main routine: find annihilation zeros, refine them, and classify true zeta zeros."""
    primes = generate_primes(num_primes)
    seeds = find_imaginary_zeros(gamma_range, primes, resolution)

    refined = []
    idx = 1

    with open(outfile, 'w') as data:
        for gamma in seeds:
            r = refine_zeta_zero(gamma)
            if r:
                z = zeta(mpc(mpf('0.5'), r))
                if (abs(z) < tol and r > 0 and gamma_range[0] <= r <= gamma_range[1] and all(abs(r - r_prev) > mpf('1e-30') for r_prev in refined)):
                    refined.append(r)
                    msg = f"{idx}: Zeta Zero {r}"
                    print(msg)
                    data.write(msg + "\n")
                    idx += 1


if __name__ == "__main__":
    if len(argv) < 2:
        raise ValueError("Usage: python zeta-util.py <gamma_upper_bound>")
    gamma_upper_bound = int(argv[1])
    compute_zeta_zeros_with_classification(
        gamma_range=(10, gamma_upper_bound),
        num_primes=10,
        resolution=5000,
        tol=1e-10,
        outfile='zeta_zeros.txt'
    )

