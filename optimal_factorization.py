from math import prod, log10, sqrt, ceil

def find_factorizations(x, beg=2, current=None, results=None):
    """
    Recursively finds all factorizations of integer x into a list of integers >= beg.
    """
    if current is None:
        current = []
    if results is None:
        results = []
    if x == 1:
        results.append(list(reversed(current)))
        return results
    for i in range(beg, x + 1):
        if x % i == 0:
            current.append(i)
            find_factorizations(x // i, i, current, results)
            current.pop()
    return results

def D_infinity(delta_p, delta_s):
    """
    Calculates the asymptotic FIR filter order parameter D∞ from ripple values.
    """
    log_dp = log10(delta_p)
    log_ds = log10(delta_s)
    return ((5.309e-3 * log_dp**2 + 7.114e-2 * log_dp - 0.4761) * log_ds
            - (2.66e-3 * log_dp**2 + 0.5941 * log_dp + 0.4278))

def L_beta(delta_p, delta_s, df):
    """
    Estimates the length of the TMFS.
    """
    return 2 * ceil((-20 * log10(sqrt(delta_p * delta_s)) - 8.4) / (30.4 * df))

def M_beta(delta_p, delta_s, df):
    """
    Estimates the polynomial order of the TMFS.
    """
    num = -20 * log10(sqrt(delta_p * delta_s)) - 8.4
    return ceil((2 * num / (30.4 * df))**0.7 - 1)

def find_optimal_factorization(delta_p, delta_s, f_p, R, F_in):
    """
    Finds the optimal factorization of decimation ratio R that minimizes
    the overall computational load for a TMFS + FIR cascade.
    """
    ret = None
    result = float("inf")

    # Case 1: Integer-only FIR cascade
    if isinstance(R, int):
        all_factorizations = find_factorizations(R)
        for n in all_factorizations:
            K = len(n)
            tmp = 0
            product = 1
            for k in range(K):
                product *= n[k]
                tmp += D_infinity(delta_p / K, delta_s) / ((2 - product / R * (1 + f_p)) / n[k] * product)
            if tmp < result:
                result = tmp
                ret = n

    # Case 2: TMFS + FIR cascade
    for N in range(2, int(R) + 1):
        all_factorizations = find_factorizations(N)
        for n in all_factorizations:
            beta = R / N
            K = len(n)
            l = L_beta(delta_p / (K + 1), delta_s, 1 - (1 + f_p) / (2 * N))
            m = M_beta(delta_p / (K + 1), delta_s, 1 - (1 + f_p) / (2 * N))
            tmp = m + l / 2 * (m + 1) / beta
            product = 1
            for k in range(K):
                product *= n[k]
                tmp += D_infinity(delta_p / (K + 1), delta_s) / ((2 - product / N * (1 + f_p)) / n[k] * beta * product)
            if tmp < result:
                result = tmp
                ret = n

    return F_in * result, ret


if __name__ == "__main__":
    # Example usage
    delta_p = 0.01
    delta_s = 0.001
    f_p = 0.9
    F_in = 1
    R = 97

    cost, factors = find_optimal_factorization(delta_p, delta_s, f_p, R, F_in)
    print(f"Optimal factorization for R = {R}:")
    print(f"Computation cost: {cost:.4f}")
    print(f"FIR factors: {factors}")
